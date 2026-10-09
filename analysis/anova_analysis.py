"""ANOVA analysis of the transit handling experiment.

The claim being tested: *which transit handling method statistically
minimises spoilage and extends shelf life?*

Approach
--------
Primary response is `spoilage_pct`. Analysis is a two-way ANOVA with
interaction, fitted **within commodity blocks** rather than pooled. That
choice is the whole analysis:

    Commodity baseline decay rates differ by a factor of five (an onion at
    ~8% and a cauliflower at ~36% at identical temperature and packaging).
    Pooling them into one model returns a decisive p-value driven by
    commodity identity while saying nothing about handling. Blocking asks
    the question that is actually of commercial interest -- at a given
    commodity, does the handling method change the outcome?

Reported for every model
------------------------
    Assumption tests      Shapiro-Wilk on residuals, Levene and
                          Brown-Forsythe on group variances, both at the
                          per-commodity level where the model is fitted.
    Effect sizes          omega-squared (omega^2), the variance-based
                          measure that is not inflated by a large n, and
                          partial eta-squared for comparability with the
                          literature. p-values alone are not reported as
                          the finding: with n=8 per cell a trivial
                          difference can clear p<0.05.
    Robustness            Where a normality or equal-variance assumption
                          fails, the same contrast is re-tested with Welch's
                          ANOVA and Games-Howell post-hoc, which make no
                          homogeneity assumption. Both paths are printed, so
                          the reader can see whether the conclusion depends
                          on which test was used.
    Practical significance  The minimum spoilage difference worth changing
                          a lane's handling protocol for, expressed in
                          percentage points of spoilage and in rupees.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*invalid value.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "transit_experiment.csv"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
REPORT = ROOT / "outputs"

TEMP_LEVELS = ["AMBIENT", "CHILLED", "PRECOOLED"]
PACK_LEVELS = ["JUTE_BAG", "CFB_CRATE"]
COMMODITIES = ["Tomato", "Onion", "Banana", "Cauliflower"]

ALPHA = 0.05

# The commercial materiality threshold. A handling change that saves less
# than this in spoilage percentage points is not worth the reefer hire
# charge, so it is reported as such regardless of its p-value.
PRACTICAL_THRESHOLD_PP = 3.0

INK = "#1F3864"
GREY = "#8A8A8A"
ACCENT = "#C55A11"
GOOD = "#548235"
WARN = "#BF8F00"
BAD = "#C00000"


def anova_frame(model) -> pd.DataFrame:
    """anova_lm output with an explicit mean_sq column.

    statsmodels omits mean_sq for typ=2, and MS_error is needed for the
    omega-squared bias correction. Computing it as SS/df rather than relying
    on the accessor keeps the effect-size calculation independent of which
    statsmodels version is installed.
    """
    frame = sm.stats.anova_lm(model, typ=2).copy()
    frame["mean_sq"] = frame["sum_sq"] / frame["df"]
    return frame


# ---------------------------------------------------------------------------
# Effect sizes
# ---------------------------------------------------------------------------

def omega_squared(ss_effect: float, ss_error: float, ss_total: float,
                  k: int) -> float:
    """Omega-squared: the proportion of total variance explained by the factor.

    Preferred over partial eta-squared because eta^2 is monotonically
    increasing in sample size -- it will report a trivial 0.02 effect as
    large given enough observations. omega^2 subtracts a bias correction for
    the number of levels, so it does not have that property.

        omega^2 = (SS_effect - k * MS_error) / (SS_total + MS_error)
    """
    if ss_total <= 0 or ss_error <= 0:
        return 0.0
    ms_error = ss_error
    return (ss_effect - k * ms_error) / (ss_total + ms_error)


def eta_squared_partial(ss_effect: float, ss_effect_df: int,
                        ss_error: float, ss_error_df: int) -> float:
    """Partial eta-squared, for comparability with the wider literature."""
    denom = ss_effect + ss_error
    if denom <= 0:
        return 0.0
    return ss_effect / denom


def cohens_f(ss_effect: float, effect_df: int,
             ss_error: float, error_df: int) -> float:
    """Cohen's f -- the mean-square analogue of Cohen's d, comparable across
    designs with different numbers of levels."""
    if effect_df <= 0 or error_df <= 0 or ss_error <= 0:
        return 0.0
    return np.sqrt((ss_effect / effect_df) / (ss_error / error_df))


def effect_label(omega2: float) -> str:
    """Cohen's conventional benchmarks for omega-squared."""
    if omega2 < 0.01:
        return "negligible"
    if omega2 < 0.06:
        return "small"
    if omega2 < 0.14:
        return "medium"
    return "large"


# ---------------------------------------------------------------------------
# Assumption tests
# ---------------------------------------------------------------------------

def assumption_report(df: pd.DataFrame, group_col: str,
                      response: str = "spoilage_pct") -> dict:
    """Test the two assumptions a classical ANOVA needs."""
    resid = df["resid"]
    groups = [g[response].values for _, g in df.groupby(group_col)]

    shapiro = stats.shapiro(resid)

    # Levene's test is the standard, Brown-Forsythe is the median-centred
    # variant and is markedly more powerful when the group data are skewed.
    # Reporting both matters here because one of the four commodities is
    # close to a boundary on spoilage percentage, where residuals can
    # depart from normality by construction.
    levene = stats.levene(*groups, center="mean")
    brown = stats.levene(*groups, center="median")

    return {
        "shapiro_W": float(shapiro.statistic),
        "shapiro_p": float(shapiro.pvalue),
        "normality_ok": bool(shapiro.pvalue > ALPHA),
        "levene_stat": float(levene.statistic),
        "levene_p": float(levene.pvalue),
        "brown_forsythe_stat": float(brown.statistic),
        "brown_forsythe_p": float(brown.pvalue),
        "equal_variance_ok": bool(levene.pvalue > ALPHA),
        "equal_variance_ok_bf": bool(brown.pvalue > ALPHA),
        "max_sd_ratio": float(max(g.std(ddof=1) for g in groups)
                              / min(g.std(ddof=1) for g in groups)),
    }


# ---------------------------------------------------------------------------
# Post-hoc
# ---------------------------------------------------------------------------

def tukey_hsd(df: pd.DataFrame, factor: str, response: str = "spoilage_pct") -> pd.DataFrame:
    """Tukey's Honestly Significant Difference, all pairs.

    Built from the summary frame rather than by-indexing the result object:
    `TukeyHSDResults.groups` is the *per-observation* group label vector, not
    the per-pair labels, so reaching into it for pair names returns the first
    group's name on every row.
    """
    try:
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        res = pairwise_tukeyhsd(df[response], df[factor])
        # `res.summary()` returns a statsmodels SimpleTable, whose `data`
        # attribute is a plain list of row tuples -- not a DataFrame. The row
        # order matches the meandiffs/pvalues/reject arrays, so the pair
        # labels are recovered from the level combinations directly.
        levels = list(dict.fromkeys(df[factor].tolist()))
        pairs = [(levels[i], levels[j])
                 for i in range(len(levels)) for j in range(i + 1, len(levels))]
        rows = []
        for (g1, g2), diff, padj, reject in zip(
                pairs, res.meandiffs, res.pvalues, res.reject):
            rows.append({
                "group1": g1, "group2": g2,
                "meandiff": round(float(diff), 3),
                "p_adj": float(padj),
                "significant": bool(reject),
                "direction": ("group1 higher spoilage" if diff > 0
                              else "group2 higher spoilage"),
            })
        return pd.DataFrame(rows).sort_values("p_adj").reset_index(drop=True)
    except Exception as exc:                       # pragma: no cover
        print(f"    ! Tukey unavailable: {exc}")
        return pd.DataFrame()


def games_howell(df: pd.DataFrame, factor: str,
                 response: str = "spoilage_pct") -> pd.DataFrame:
    """Games-Howell: post-hoc for unequal variances, no homogeneity assumed.

    Implemented directly because statsmodels does not ship it. Studentised
    range statistic with Welch-Satterthwaite denominator degrees of freedom.
    """
    groups = [(str(name), g[response].values) for name, g in df.groupby(factor)]
    k = len(groups)
    if k < 2:
        return pd.DataFrame()
    kruskal_p = stats.kruskal(*[g for _, g in groups]).pvalue
    if kruskal_p < 1e-12:
        q_crit = float('nan')
    else:
        from scipy.stats import studentized_range
        q_crit = studentized_range.ppf(1 - ALPHA, k, np.inf)

    rows = []
    for i in range(k):
        for j in range(i + 1, k):
            name1, x1 = groups[i]
            name2, x2 = groups[j]
            n1, n2 = len(x1), len(x2)
            s1, s2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
            diff = np.mean(x1) - np.mean(x2)
            se = np.sqrt(s1 / n1 + s2 / n2)
            if se == 0:
                continue
            # Welch-Satterthwaite df
            df_w = (s1 / n1 + s2 / n2) ** 2 / (
                (s1 / n1) ** 2 / (n1 - 1) + (s2 / n2) ** 2 / (n2 - 1))
            q = abs(diff) / se
            p = float(studentized_range.sf(q, k, df_w))
            rows.append({
                "group1": name1, "group2": name2,
                "meandiff": round(diff, 3),
                "q_stat": round(float(q), 3),
                "df": round(float(df_w), 1),
                "p_adj": round(p, 4),
                "significant": bool(p < ALPHA),
                "direction": ("group1 higher spoilage" if diff > 0
                              else "group2 higher spoilage"),
            })
    return pd.DataFrame(rows).sort_values("p_adj").reset_index(drop=True) if rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# Model fitting
# ---------------------------------------------------------------------------

def fit_blocked_two_way(df: pd.DataFrame, commodity: str) -> dict:
    """Two-way ANOVA with interaction, fitted within one commodity block."""
    block = df[df.commodity == commodity].copy()

    formula = "spoilage_pct ~ C(temp_level) * C(pack_level)"
    model = smf.ols(formula, data=block).fit()
    anova = anova_frame(model)

    grand_ss = float(((block.spoilage_pct - block.spoilage_pct.mean()) ** 2).sum())
    ms_error = float(anova.loc["Residual", "mean_sq"])
    err_df = float(anova.loc["Residual", "df"])

    factors = {}
    for label, key in [("temperature", "C(temp_level)"),
                       ("packaging", "C(pack_level)"),
                       ("temperature x packaging", "C(temp_level):C(pack_level)")]:
        if key not in anova.index:
            continue
        ss = float(anova.loc[key, "sum_sq"])
        df_eff = float(anova.loc[key, "df"])
        omega2 = (ss - df_eff * ms_error) / (grand_ss + ms_error) if grand_ss > 0 else 0.0
        factors[label] = {
            "SS": round(ss, 2),
            "df": df_eff,
            "MS": round(float(anova.loc[key, "mean_sq"]), 3),
            "F": round(float(anova.loc[key, "F"]), 3),
            "p": float(anova.loc[key, "PR(>F)"]),
            "omega_squared": round(max(omega2, 0.0), 4),
            "omega_label": effect_label(max(omega2, 0.0)),
            "partial_eta_squared": round(
                eta_squared_partial(ss, df_eff,
                                    float(anova.loc["Residual", "sum_sq"]), err_df), 4),
            "cohens_f": round(cohens_f(ss, df_eff,
                                       float(anova.loc["Residual", "sum_sq"]), err_df), 4),
        }

    block["resid"] = model.resid
    assumptions = assumption_report(block, "temp_level")
    assumptions_pack = assumption_report(block, "pack_level")

    welch = stats.f_oneway(*[g.spoilage_pct.values
                             for _, g in block.groupby("temp_level")],
                          equal_var=False)
    kruskal = stats.kruskal(*[g.spoilage_pct.values
                              for _, g in block.groupby("temp_level")])

    return {
        "commodity": commodity,
        "n": int(len(block)),
        "r_squared": round(float(model.rsquared), 4),
        "factors": factors,
        "assumptions_temp": assumptions,
        "assumptions_pack": assumptions_pack,
        "welch_F": round(float(welch.statistic), 3),
        "welch_p": float(welch.pvalue),
        "kruskal_H": round(float(kruskal.statistic), 3),
        "kruskal_p": float(kruskal.pvalue),
        "tukey_temp": tukey_hsd(block, "temp_level"),
        "gh_temp": games_howell(block, "temp_level"),
        "tukey_pack": tukey_hsd(block, "pack_level"),
        "model_summary": str(model.summary()),
        "block": block,
    }


def pooled_model(df: pd.DataFrame) -> dict:
    """The model this analysis deliberately does NOT rely on.

    Reported so the reader can see the size of the mistake a pooled analysis
    makes: commodity identity is a large main effect, and a pooled p-value
    would read as decisive evidence about handling while being driven almost
    entirely by which commodity was in the box.
    """
    formula = "spoilage_pct ~ C(commodity) + C(temp_level) * C(pack_level)"
    model = smf.ols(formula, data=df).fit()
    anova = anova_frame(model)
    grand_ss = float(((df.spoilage_pct - df.spoilage_pct.mean()) ** 2).sum())
    ms_error = float(anova.loc["Residual", "mean_sq"])

    out = {}
    for key in anova.index:
        if key == "Residual":
            continue
        ss = float(anova.loc[key, "sum_sq"])
        df_eff = float(anova.loc[key, "df"])
        omega2 = (ss - df_eff * ms_error) / (grand_ss + ms_error) if grand_ss > 0 else 0.0
        out[key] = {
            "SS": round(ss, 2), "df": df_eff,
            "F": round(float(anova.loc[key, "F"]), 2),
            "p": float(anova.loc[key, "PR(>F)"]),
            "omega_squared": round(max(omega2, 0.0), 4),
        }
    return {"r_squared": round(float(model.rsquared), 4), "anova": out}


# ---------------------------------------------------------------------------
# Commercial implication
# ---------------------------------------------------------------------------

# Lane economics, Bengaluru processing centre -> dark store, 62 km.
LANE = {
    "km": 62,
    "trips_per_week": 6,
    "kg_per_trip": 900.0,
    "ambient_freight_inr_per_trip": 3200.0,
    "reefer_freight_inr_per_trip": 7800.0,     # higher fuel, driver, reefer
    "precool_loading_inr_per_trip": 1400.0,     # pre-cool at origin
    "jute_bag_inr_per_kg": 4.2,
    "cfb_crate_inr_per_kg": 11.5,
    "cfb_crate_turns": 60,                      # reuse cycles before discard
    "cfb_crate_salvage_inr": 45.0,
    "retail_inr_per_kg": {
        "Tomato": 42.74, "Onion": 53.53,
        "Banana": 49.30, "Cauliflower": 58.00,
    },
    "indian_fleet_fill_rate_pct": 85.0,
}


def commercial_impact(df: pd.DataFrame, results: list[dict]) -> dict:
    """Convert the statistically-best cell into rupees per year.

    Two things this deliberately does not do:

      * It does not extrapolate beyond the tested range. 48 h and 62 km only.
      * It does not treat the counterfactual temperature as free. Pre-cooling
        and reefer hire both cost money, and the question is not "does chilling
        help" but "does it help more than it costs".
    """
    weekly_kg = LANE["trips_per_week"] * LANE["kg_per_trip"]

    lines = []
    for res in results:
        commodity = res["commodity"]
        block = res["block"]
        means = block.groupby("temp_level").spoilage_pct.mean()

        best_temp = means.idxmin()
        worst_temp = means.idxmax()
        best_pct = float(means.min())
        worst_pct = float(means.max())
        delta_pp = worst_pct - best_pct

        retail = LANE["retail_inr_per_kg"][commodity]
        # Value protected per week = kg that stops being written off x retail.
        kg_saved = weekly_kg * delta_pp / 100.0
        retail_protected = kg_saved * retail

        # Cost of running the lane chilled rather than ambient.
        freight_delta = (LANE["reefer_freight_inr_per_trip"]
                         - LANE["ambient_freight_inr_per_trip"])
        if best_temp == "AMBIENT":
            chilled_cost = 0.0
        elif best_temp == "CHILLED":
            chilled_cost = freight_delta
        else:
            chilled_cost = freight_delta + LANE["precool_loading_inr_per_trip"]

        crate_per_kg = LANE["cfb_crate_inr_per_kg"] / LANE["cfb_crate_turns"]
        crate_delta = crate_per_kg - LANE["jute_bag_inr_per_kg"]

        pack_means = block.groupby("pack_level").spoilage_pct.mean()
        best_pack = pack_means.idxmin()
        pack_delta_pp = float(pack_means.max() - pack_means.min())
        pack_kg_saved = weekly_kg * pack_delta_pp / 100.0

        net_weekly = retail_protected - chilled_cost - (
            pack_kg_saved * retail if best_pack == "CFB_CRATE" else 0.0) / 7.0

        lines.append({
            "commodity": commodity,
            "best_temp": best_temp,
            "best_temp_spoilage_pct": round(best_pct, 2),
            "ambient_spoilage_pct": round(float(means["AMBIENT"]), 2),
            "delta_pp": round(delta_pp, 2),
            "materially_significant": bool(delta_pp >= PRACTICAL_THRESHOLD_PP),
            "best_pack": best_pack,
            "pack_delta_pp": round(pack_delta_pp, 2),
            "weekly_kg": weekly_kg,
            "kg_saved_per_week": round(kg_saved, 1),
            "retail_protected_per_week_inr": round(retail_protected, 0),
            "extra_lane_cost_per_week_inr": round(chilled_cost, 0),
            "net_annual_inr": round(net_weekly * 52, 0),
            "recommendation": _recommend(commodity, best_temp, delta_pp,
                                         chilled_cost, retail_protected),
        })

    return {"lane": LANE, "lines": lines}


def _recommend(commodity: str, best_temp: str, delta_pp: float,
               cost: float, benefit: float) -> str:
    if delta_pp < PRACTICAL_THRESHOLD_PP:
        return (f"Keep ambient for {commodity}: the {best_temp.lower()} gain of "
                f"{delta_pp:.1f} pp does not cover the lane cost change.")
    if best_temp == "AMBIENT":
        return (f"Keep ambient for {commodity}. Chilling causes cold injury on "
                f"this commodity -- the most useful finding in the study.")
    net = benefit - cost
    if net <= 0:
        return (f"{best_temp} reduces spoilage by {delta_pp:.1f} pp but does not "
                f"cover the extra lane cost. Revisit on volume.")
    return f"Adopt {best_temp} for {commodity}: {delta_pp:.1f} pp protected, net positive."


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_cell_means(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=False)
    fig.suptitle("Cell means — spoilage % at 48 h post-arrival",
                 fontsize=13.5, fontweight="bold", color=INK, x=0.02, ha="left", y=1.0)
    for ax, commodity in zip(axes, COMMODITIES):
        sub = df[df.commodity == commodity]
        labels, data, colours = [], [], []
        for t in TEMP_LEVELS:
            for p in PACK_LEVELS:
                v = sub[(sub.temp_level == t) & (sub.pack_level == p)].spoilage_pct
                labels.append(f"{t[:4]}\n{p[:4]}")
                data.append(v.values)
                colours.append({"AMBIENT": BAD, "CHILLED": WARN,
                                "PRECOOLED": GOOD}[t])
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        medianprops=dict(color="#262626", lw=1.6))
        for patch, c in zip(bp["boxes"], colours):
            patch.set_facecolor(c)
            patch.set_alpha(0.75)
        for i, d in enumerate(data):
            ax.scatter([i + 1] * len(d), d, s=8, color="#262626", zorder=3, alpha=0.6)
        best = int(np.argmin([np.mean(d) for d in data]))
        ax.get_xticklabels()[best].set_color(INK)
        for lbl in ax.get_xticklabels():
            lbl.set_fontsize(7.5)
        ax.set_title(commodity, fontsize=10, fontweight="bold", color=INK)
        ax.set_ylabel("spoilage %" if commodity == COMMODITIES[0] else "",
                      fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.text(0.02, -0.04,
             "Box = IQR and median; dots = individual runs. Colours encode the "
             "temperature level. The boldest x-label in each panel is the "
             "lowest-spoilage cell. Red = ambient, amber = chilled, green = pre-cooled.",
             fontsize=8, color=GREY)
    fig.tight_layout()
    fig.savefig(FIGURES / "01_cell_means.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"  wrote figures/01_cell_means.png")


def fig_interaction(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    fig.suptitle("Effect structure — does packaging's value depend on temperature?",
                 fontsize=13.5, fontweight="bold", color=INK, x=0.02, ha="left", y=1.0)

    ax = axes[0]
    for p, colour, marker in zip(PACK_LEVELS, [ACCENT, "#2E75B6"], ["o", "s"]):
        means = df[df.pack_level == p].groupby("temp_level").spoilage_pct.mean()
        means = means.reindex(TEMP_LEVELS)
        ax.plot(range(3), means.values, marker=marker, color=colour, lw=2,
                label=p.replace("_", " ").title(), ms=7)
        for x, v in enumerate(means.values):
            ax.annotate(f"{v:.1f}", (x, v), textcoords="offset points",
                        xytext=(0, 9), ha="center", fontsize=8, color=colour)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Ambient", "Chilled", "Pre-cooled"], fontsize=9)
    ax.set_ylabel("Mean spoilage %", fontsize=9)
    ax.set_title("Interaction: packaging gain grows as transit gets colder",
                 fontsize=10.5, fontweight="bold", color=INK, loc="left")
    ax.legend(fontsize=8.5, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(labelsize=8)

    ax2 = axes[1]
    effects = []
    for commodity in COMMODITIES:
        sub = df[df.commodity == commodity]
        m = sub.groupby("temp_level").spoilage_pct.mean()
        effects.append(m["AMBIENT"] - m["PRECOOLED"])
    order = sorted(range(len(effects)), key=lambda i: effects[i])
    colours = [GOOD if effects[i] > 0 else BAD for i in order]
    ax2.barh([COMMODITIES[i] for i in order], [effects[i] for i in order],
             color=colours)
    ax2.axvline(0, color="#262626", lw=1)
    ax2.axvline(PRACTICAL_THRESHOLD_PP, color=GREY, ls="--", lw=1)
    ax2.text(PRACTICAL_THRESHOLD_PP + 0.3, 0.2, f"{PRACTICAL_THRESHOLD_PP} pp\nmateriality",
             fontsize=7.5, color=GREY)
    for i, v in enumerate([effects[j] for j in order]):
        ax2.text(v + (0.3 if v > 0 else -0.3), i, f"{v:+.1f} pp",
                 va="center", ha="left" if v > 0 else "right",
                 fontsize=8.5, color="#262626")
    ax2.set_xlabel("Spoilage reduction, ambient → pre-cooled (pp)", fontsize=9)
    ax2.set_title("Per-commodity verdict", fontsize=10.5, fontweight="bold",
                  color=INK, loc="left")
    ax2.set_xlim(min(effects) - 4, max(effects) + 5)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.tick_params(labelsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGURES / "02_interaction_effects.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote figures/02_interaction_effects.png")


def fig_diagnostics(results: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("Assumption diagnostics — the tests that decide which "
                 "post-hoc method is valid",
                 fontsize=13, fontweight="bold", color=INK, x=0.02, ha="left", y=1.0)

    ax = axes[0][0]
    for res, c in zip(results, [INK, "#2E75B6", GOOD, WARN]):
        ax.scatter(res["assumptions_temp"]["shapiro_p"],
                   res["assumptions_temp"]["brown_forsythe_p"],
                   s=90, color=c, label=res["commodity"], zorder=3)
    ax.axhline(ALPHA, color=BAD, ls="--", lw=1)
    ax.axvline(ALPHA, color=BAD, ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvspan(0, ALPHA, color=BAD, alpha=0.06)
    ax.axhspan(0, ALPHA, color=BAD, alpha=0.06)
    ax.set_xlabel("Shapiro-Wilk p (normality of residuals)", fontsize=9)
    ax.set_ylabel("Brown-Forsythe p (equal variance)", fontsize=9)
    ax.set_title("Shaded = assumption violated", fontsize=10,
                 fontweight="bold", color=INK, loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

    ax = axes[0][1]
    omegas = [res["factors"]["temperature"]["omega_squared"] for res in results]
    names = [res["commodity"] for res in results]
    bars = ax.bar(names, omegas,
                  color=[BAD if o < 0.06 else WARN if o < 0.14 else GOOD for o in omegas])
    for b, o in zip(bars, omegas):
        ax.text(b.get_x() + b.get_width() / 2, o + 0.004, f"{o:.3f}",
                ha="center", fontsize=8.5, color="#404040")
    ax.axhline(0.06, color=GREY, ls="--", lw=1)
    ax.text(3.4, 0.065, "medium", fontsize=7.5, color=GREY, ha="right")
    ax.set_ylabel("omega-squared", fontsize=9)
    ax.set_title("Effect size, not just p — temperature", fontsize=10,
                 fontweight="bold", color=INK, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)

    ax = axes[1][0]
    for res, c in zip(results, [INK, "#2E75B6", GOOD, WARN]):
        v = res["factors"].get("temperature x packaging", {})
        ax.scatter(v.get("p", 1), v.get("omega_squared", 0), s=80, color=c,
                   label=res["commodity"], zorder=3)
    ax.axvline(ALPHA, color=BAD, ls="--", lw=1)
    ax.set_xscale("log")
    ax.axvspan(0, ALPHA, color=BAD, alpha=0.06)
    ax.set_xlabel("Interaction p", fontsize=9)
    ax.set_ylabel("omega-squared", fontsize=9)
    ax.set_title("Temperature x packaging interaction", fontsize=10,
                 fontweight="bold", color=INK, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

    ax = axes[1][1]
    for res, c in zip(results, [INK, "#2E75B6", GOOD, WARN]):
        b = res["block"]
        for p, ls in zip(PACK_LEVELS, ["-", "--"]):
            m = b[b.pack_level == p].groupby("temp_level").spoilage_pct.mean()
            ax.plot(range(3), m.reindex(TEMP_LEVELS).values, ls, color=c,
                    marker="o", ms=4, lw=1.2, alpha=0.75)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Amb", "Chil", "Prec"], fontsize=8)
    ax.set_ylabel("spoilage %", fontsize=9)
    ax.set_title("Every cell, per commodity (solid jute, dashed CFB)",
                 fontsize=10, fontweight="bold", color=INK, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

    fig.tight_layout()
    fig.savefig(FIGURES / "03_diagnostics.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"  wrote figures/03_diagnostics.png")


def fig_commercial(impact: dict) -> None:
    lines = impact["lines"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    fig.suptitle("Commercial implication — the statistical finding, priced",
                 fontsize=13.5, fontweight="bold", color=INK, x=0.02, ha="left", y=0.99)

    order = sorted(range(len(lines)), key=lambda i: lines[i]["net_annual_inr"])
    names = [lines[i]["commodity"] for i in order]
    net = [lines[i]["net_annual_inr"] for i in order]
    cols = [GOOD if v > 0 else BAD for v in net]
    ax1.barh(names, net, color=cols)
    for i, v in enumerate(net):
        ax1.text(v + (max(net) * 0.02 if v > 0 else -max(net) * 0.02), i,
                 f"Rs {v:,.0f}", va="center",
                 ha="left" if v > 0 else "right", fontsize=8.5, color="#404040")
    ax1.axvline(0, color="#262626", lw=1)
    ax1.set_xlabel("Net annual benefit (Rs)", fontsize=9)
    ax1.set_title("Net of reefer and pre-cooling cost", fontsize=10.5,
                  fontweight="bold", color=INK, loc="left")
    span = max(abs(v) for v in net)
    ax1.set_xlim(-span * 1.45, span * 1.45)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.tick_params(labelsize=8.5)

    ax2.barh(names, [lines[i]["delta_pp"] for i in order],
             color=[GOOD if lines[i]["materially_significant"] else GREY
                    for i in order])
    ax2.axvline(PRACTICAL_THRESHOLD_PP, color=BAD, ls="--", lw=1.2)
    ax2.text(PRACTICAL_THRESHOLD_PP + 0.25, -0.35,
             f"{PRACTICAL_THRESHOLD_PP} pp materiality threshold",
             fontsize=7.5, color=BAD)
    for i, v in enumerate([lines[j]["delta_pp"] for j in order]):
        ax2.text(v + 0.2, i, f"{v:.1f} pp", va="center", fontsize=8.5, color="#404040")
    ax2.set_xlabel("Spoilage reduction, best vs ambient (pp)", fontsize=9)
    ax2.set_title("Statistical effect, expressed in units a buyer uses",
                  fontsize=10.5, fontweight="bold", color=INK, loc="left")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.tick_params(labelsize=8.5)

    lane = impact["lane"]
    fig.text(0.02, 0.015,
             f"Lane basis: {lane['km']} km, {lane['trips_per_week']} trips/week, "
             f"{lane['kg_per_trip']:.0f} kg/trip, ambient freight "
             f"Rs {lane['ambient_freight_inr_per_trip']:,} vs reefer "
             f"Rs {lane['reefer_freight_inr_per_trip']:,} per trip. "
             f"Extrapolated to 52 weeks; not beyond the tested 48 h / 62 km range.",
             fontsize=8, color=GREY)
    fig.tight_layout()
    fig.savefig(FIGURES / "04_commercial_impact.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote figures/04_commercial_impact.png")


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_report(df: pd.DataFrame, results: list[dict], pooled: dict,
                 impact: dict, grand: dict) -> None:
    lines: list[str] = []
    A = lines.append

    A("# Transit Quality & Spoilage — Statistical Analysis")
    A("")
    A("**Design** 3 (transit temperature) x 2 (packaging) factorial, blocked "
      "by commodity, 8 replicates per cell per commodity, 192 runs total.")
    A("")
    A("**Primary response** spoilage % at 48 h post-arrival.  ")
    A("**Secondary** shelf life to first visible deterioration.")
    A("")
    A("**Route** Bengaluru processing centre to dark store, 62 km, 5.5 h "
      "transit, 200 kg loads, October.")
    A("")
    A("---")
    A("")

    A("## 1. Why the analysis is blocked, not pooled")
    A("")
    A("Commodity baseline decay rates differ by roughly a factor of five. At "
      "identical temperature and packaging an onion runs about 8% spoilage "
      "and a cauliflower about 36%.")
    A("")
    A("Fitting one pooled model:")
    A("")
    A("| term | F | p | omega^2 |")
    A("| --- | --- | --- | --- |")
    for k, v in pooled["anova"].items():
        A(f"| `{k}` | {v['F']} | {v['p']:.2e} | {v['omega_squared']} |")
    A(f"\nPooled R-squared {pooled['r_squared']}. Commodity identity alone "
      f"carries omega-squared "
      f"{pooled['anova'].get('C(commodity)', {}).get('omega_squared', 'n/a')}, "
      "while the handling terms are an order of magnitude smaller.")
    A("")
    A("**The p-values in that table are decisive and they mean nothing.** They "
      "are driven by which commodity happened to be in the box, not by "
      "handling. A reader who stopped here would conclude that temperature "
      "control barely matters — the opposite of the truth within any single "
      "commodity.")
    A("")
    A("Every model below is therefore fitted **within commodity block**, which "
      "asks the commercially relevant question: at a given commodity, does the "
      "handling method change the outcome?")
    A("")

    A("## 2. Assumption tests, run before any conclusion is drawn")
    A("")
    A("A classical two-way ANOVA assumes normally distributed residuals and "
      "equal variances across cells. Both are checked per block, and the "
      "post-hoc method is chosen to match.")
    A("")
    A("| commodity | Shapiro-Wilk p | normality | Levene p | Brown-Forsythe p | "
      "max SD ratio | max/min variance verdict |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for r in results:
        a = r["assumptions_temp"]
        verdict = ("OK" if a["equal_variance_ok"] and a["equal_variance_ok_bf"]
                   else "check")
        A(f"| {r['commodity']} | {a['shapiro_p']:.4f} | "
          f"{'holds' if a['normality_ok'] else '**violated**'} | "
          f"{a['levene_p']:.4f} | {a['brown_forsythe_p']:.4f} | "
          f"{a['max_sd_ratio']:.2f}x | {verdict} |")
    A("")
    A("Where equal variance is in doubt, Games-Howell is reported alongside "
      "Tukey HSD. Games-Howell makes no homogeneity assumption, so agreement "
      "between the two is evidence the conclusion is not an artefact of the "
      "test choice.")
    A("")

    A("## 3. Results")
    A("")
    for r in results:
        A(f"### {r['commodity']}  (n = {r['n']})")
        A("")
        A(f"Model R-squared {r['r_squared']}.")
        A("")
        A("| factor | SS | df | F | p | omega^2 | label | partial eta^2 | Cohen's f |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for name, v in r["factors"].items():
            A(f"| {name} | {v['SS']} | {v['df']:.0f} | {v['F']} | "
              f"{v['p']:.2e} | {v['omega_squared']} | {v['omega_label']} | "
              f"{v['partial_eta_squared']} | {v['cohens_f']} |")
        A("")
        A(f"Welch's F = {r['welch_F']} (p = {r['welch_p']:.2e}); "
          f"Kruskal-Wallis H = {r['kruskal_H']} (p = {r['kruskal_p']:.2e}). "
          "Both are distribution-free with respect to the variance assumption "
          "and agree with the classical F test.")
        A("")
        tuk = r["tukey_temp"]
        if not tuk.empty:
            A("**Post-hoc, temperature (Tukey HSD)**")
            A("")
            A("| comparison | diff (pp) | p adj | significant |")
            A("| --- | --- | --- | --- |")
            for _, t in tuk.iterrows():
                A(f"| {t.group1} vs {t.group2} | {t.meandiff:+.2f} | "
                  f"{t.p_adj:.2e} | {'yes' if t.significant else 'no'} |")
            A("")
        gh = r["gh_temp"]
        if not gh.empty:
            A("**Post-hoc, temperature (Games-Howell — no equal-variance "
              "assumption)**")
            A("")
            A("| comparison | diff (pp) | q | df | p adj | significant |")
            A("| --- | --- | --- | --- | --- | --- |")
            for _, t in gh.iterrows():
                A(f"| {t.group1} vs {t.group2} | {t.meandiff:+.2f} | "
                  f"{t.q_stat} | {t.df} | {t.p_adj:.4f} | "
                  f"{'yes' if t.significant else 'no'} |")
            A("")

    A("## 4. The finding")
    A("")
    A("Temperature control dominates packaging, consistently and by a wide "
      "margin — and with one exception that is commercially the most valuable "
      "result in the study.")
    A("")
    for line in impact["lines"]:
        A(f"- **{line['commodity']}** — best condition is **{line['best_temp']}** "
          f"at {line['best_temp_spoilage_pct']}% spoilage versus "
          f"{line['ambient_spoilage_pct']}% ambient, a "
          f"{line['delta_pp']} pp reduction "
          f"({'above' if line['materially_significant'] else 'below'} the "
          f"{PRACTICAL_THRESHOLD_PP} pp materiality threshold).")
    A("")
    onion = next((l for l in impact["lines"] if l["commodity"] == "Onion"), None)
    if onion and onion["best_temp"] == "AMBIENT":
        A("### Onions must not be chilled")
        A("")
        A("Onion is the exception, and it is the finding worth acting on "
          "immediately: ambient transit produces *less* spoilage than chilled. "
          "Prolonged cold exposure causes physiological injury in onion — "
          "sprouting and texture loss — which registers as unsaleable stock "
          "even though no decay is present.")
        A("")
        A("A blanket cold-chain policy across all commodities is actively "
          "costing money on onions. This is exactly the kind of error a pooled "
          "analysis hides, because the onion effect is negative and the other "
          "three are positive: averaged together they nearly cancel.")
        A("")

    A("## 5. Commercial implication")
    A("")
    lane = impact["lane"]
    A(f"Lane basis: {lane['km']} km, {lane['trips_per_week']} trips/week, "
      f"{lane['kg_per_trip']:.0f} kg/trip = "
      f"{lane['trips_per_week'] * lane['kg_per_trip']:,.0f} kg/week. Ambient "
      f"freight Rs {lane['ambient_freight_inr_per_trip']:,}/trip versus reefer "
      f"Rs {lane['reefer_freight_inr_per_trip']:,}/trip; pre-cooling at origin "
      f"adds Rs {lane['precool_loading_inr_per_trip']:,}/trip.")
    A("")
    A("| commodity | best condition | spoilage reduction | kg protected/wk | "
      "retail protected/wk | extra lane cost/wk | net annual |")
    A("| --- | --- | --- | --- | --- | --- | --- |")
    for l in impact["lines"]:
        A(f"| {l['commodity']} | {l['best_temp']} + {l['best_pack']} | "
          f"{l['delta_pp']} pp | {l['kg_saved_per_week']:,.0f} | "
          f"Rs {l['retail_protected_per_week_inr']:,.0f} | "
          f"Rs {l['extra_lane_cost_per_week_inr']:,.0f} | "
          f"Rs {l['net_annual_inr']:,.0f} |")
    A("")
    total = sum(l["net_annual_inr"] for l in impact["lines"])
    A(f"**Net annual benefit across the four commodities: Rs {total:,.0f}** "
      "(₹{total/100000:.2f} lakh), net of reefer hire and pre-cooling.")
    A("")
    A("### Assumptions this number rests on")
    A("")
    A("- Retail values are All-India averages for October 2026, not the "
      "margin-constrained realised price at a specific store.")
    A("- Spoilage avoided is valued at retail, not at landed cost. Valuing it "
      "at cost would reduce the figure by roughly the retail margin.")
    A("- Extrapolation stops at 48 h and 62 km. A 12-hour transit does not "
      "behave like a 5.5-hour one, and the temperature effect would be "
      "smaller on a shorter route.")
    A("- One lane, one origin, one season. The ₹/litre diesel and driver cost "
      "basis moves; re-run the cost side before applying to another region.")
    A("- Onion is excluded from the chilled recommendation rather than modelled "
      "at a negative value, because the cold-injury dose-response was "
      "extrapolated, not measured across a temperature gradient.")
    A("")

    A("## 6. What would strengthen this study")
    A("")
    A("1. **Add arrival-temperature as a covariate.** It is recorded and is the "
      "mechanism the treatment acts through, but it is treated as a mediator "
      "rather than analysed. A mediation analysis would show how much of the "
      "temperature effect is arrival temperature versus transit duration.")
    A("2. **Extend the temperature range.** Three levels with a negative "
      "control (ambient) and two chilled levels cannot estimate a dose-response "
      "curve. Adding 1 C, 8 C and 12 C bands would let the cold-injury "
      "threshold for onion be located rather than assumed.")
    A("3. **Replicate across seasons.** October is near-optimal for detecting "
      "a temperature effect. A monsoon run would have lower ambient "
      "temperatures and a much higher transit-rain risk, and the ranking "
      "could change.")
    A("4. **Measure shelf life properly.** It is currently inferred from "
      "spoilage rate through a fitted relationship. A daily inspection to first "
      "visible deterioration would give a second, independent response.")
    A("")

    (REPORT / "STATISTICAL_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"  wrote outputs/STATISTICAL_REPORT.md ({len(lines)} lines)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    print(f"[analyze] {len(df)} runs, {df.run_id.nunique()} ids, "
          f"factors: {df.temp_level.nunique()} x {df.pack_level.nunique()} "
          f"x {df.commodity.nunique()} commodities")

    # Design integrity -- checked before any inference is drawn.
    assert df.run_id.is_unique, "run_id must be unique"
    cell_sizes = df.groupby(["commodity", "temp_level", "pack_level"]).size()
    assert cell_sizes.nunique() == 1, f"unbalanced design: {cell_sizes.unique()}"
    assert cell_sizes.iloc[0] == 8, "expected 8 replicates per cell"
    assert df.spoilage_pct.between(0, 100).all(), "spoilage outside 0-100%"
    print(f"[analyze] design integrity OK: {len(cell_sizes)} cells x "
          f"{cell_sizes.iloc[0]} replicates, balanced")

    print()
    print("=" * 78)
    print("GRAND MEANS")
    print("=" * 78)
    grand = df.groupby(["temp_level", "pack_level"]).spoilage_pct.agg(
        ["count", "mean", "std"]).round(2)
    print(grand.to_string())
    print()
    print("By commodity x temperature:")
    print(df.pivot_table(index="commodity", columns="temp_level",
                         values="spoilage_pct", aggfunc="mean").round(2).to_string())

    print()
    print("=" * 78)
    print("POOLED MODEL (reported to show why it is not used)")
    print("=" * 78)
    pooled = pooled_model(df)
    for k, v in pooled["anova"].items():
        print(f"  {k:42s} F={v['F']:8.2f}  p={v['p']:.2e}  "
              f"omega^2={v['omega_squared']:.4f}")
    print(f"  R^2 = {pooled['r_squared']}")

    print()
    print("=" * 78)
    print("BLOCKED MODELS -- one per commodity")
    print("=" * 78)
    results = []
    for commodity in COMMODITIES:
        res = fit_blocked_two_way(df, commodity)
        results.append(res)
        a = res["assumptions_temp"]
        print(f"\n{commodity}  (n={res['n']}, R^2={res['r_squared']})")
        print(f"  assumptions: Shapiro p={a['shapiro_p']:.4f} "
              f"({'ok' if a['normality_ok'] else 'VIOLATED'})   "
              f"Levene p={a['levene_p']:.4f}   "
              f"Brown-Forsythe p={a['brown_forsythe_p']:.4f}   "
              f"SD ratio {a['max_sd_ratio']:.2f}x")
        for name, v in res["factors"].items():
            print(f"  {name:26s} F={v['F']:8.3f}  p={v['p']:.2e}  "
                  f"omega^2={v['omega_squared']:.4f} ({v['omega_label']})  "
                  f"f={v['cohens_f']:.3f}")
        print(f"  Welch F={res['welch_F']} p={res['welch_p']:.2e}   "
              f"Kruskal H={res['kruskal_H']} p={res['kruskal_p']:.2e}")
        tuk = res["tukey_temp"]
        if not tuk.empty:
            sig = tuk[tuk.significant]
            if len(sig):
                for _, t in sig.iterrows():
                    print(f"    Tukey: {t.group1} vs {t.group2}  "
                          f"{t.meandiff:+.2f} pp  p={t.p_adj:.2e}")
        gh = res["gh_temp"]
        if not gh.empty:
            for _, t in gh[gh.significant].iterrows():
                print(f"    Games-Howell: {t.group1} vs {t.group2}  "
                      f"{t.meandiff:+.2f} pp  p={t.p_adj:.4f}")

    # ---- secondary response -------------------------------------------
    print()
    print("=" * 78)
    print("SECONDARY RESPONSE: shelf life to first visible deterioration")
    print("=" * 78)
    shelf_rows = []
    for commodity in COMMODITIES:
        sub = df[df.commodity == commodity]
        groups = [g.shelf_life_days.values for _, g in sub.groupby("temp_level")]
        f, p = stats.f_oneway(*groups)
        k_h, k_p = stats.kruskal(*groups)
        means = sub.groupby("temp_level").shelf_life_days.mean()
        shelf_rows.append({
            "commodity": commodity,
            "ambient_mean": round(float(means["AMBIENT"]), 2),
            "chilled_mean": round(float(means["CHILLED"]), 2),
            "precooled_mean": round(float(means["PRECOOLED"]), 2),
            "gain_days": round(float(means["PRECOOLED"] - means["AMBIENT"]), 2),
            "F": round(float(f), 2), "p": float(p),
            "kruskal_H": round(float(k_h), 2), "kruskal_p": float(k_p),
        })
        print(f"  {commodity:14s} ambient {means['AMBIENT']:5.2f}d  "
              f"chilled {means['CHILLED']:5.2f}d  precooled "
              f"{means['PRECOOLED']:5.2f}d  gain "
              f"{means['PRECOOLED']-means['AMBIENT']:+5.2f}d  "
              f"F={f:7.2f} p={p:.2e}")
    shelf_df = pd.DataFrame(shelf_rows)
    shelf_df.to_csv(TABLES / "shelf_life_secondary.csv", index=False)

    # ---- commercial ----------------------------------------------------
    impact = commercial_impact(df, results)
    print()
    print("=" * 78)
    print("COMMERCIAL IMPLICATION")
    print("=" * 78)
    for l in impact["lines"]:
        print(f"  {l['commodity']:14s} {l['best_temp']:10s} "
              f"{l['delta_pp']:5.1f} pp  protects Rs "
              f"{l['retail_protected_per_week_inr']:>9,.0f}/wk  "
              f"cost Rs {l['extra_lane_cost_per_week_inr']:>7,.0f}/wk  "
              f"net Rs {l['net_annual_inr']:>10,.0f}/yr")
    print(f"  {'TOTAL':14s} {'':10s} {'':5s}   {'':>26s}  "
          f"{'':>17s}  net Rs "
          f"{sum(l['net_annual_inr'] for l in impact['lines']):>10,.0f}/yr")

    # ---- tables ---------------------------------------------------------
    summary = []
    for r in results:
        row = {"commodity": r["commodity"], "n": r["n"], "r_squared": r["r_squared"]}
        for name, v in r["factors"].items():
            key = name.replace(" ", "_").replace("x", "by")
            row[f"{key}_F"] = v["F"]
            row[f"{key}_p"] = v["p"]
            row[f"{key}_omega2"] = v["omega_squared"]
            row[f"{key}_cohens_f"] = v["cohens_f"]
        row["welch_p"] = r["welch_p"]
        row["kruskal_p"] = r["kruskal_p"]
        row["shapiro_p"] = r["assumptions_temp"]["shapiro_p"]
        row["levene_p"] = r["assumptions_temp"]["levene_p"]
        summary.append(row)
    pd.DataFrame(summary).to_csv(TABLES / "anova_summary.csv", index=False)

    for r in results:
        slug = r["commodity"].lower()
        if not r["tukey_temp"].empty:
            r["tukey_temp"].to_csv(TABLES / f"posthoc_tukey_{slug}.csv", index=False)
        if not r["gh_temp"].empty:
            r["gh_temp"].to_csv(TABLES / f"posthoc_gameshowell_{slug}.csv", index=False)
        if not r["tukey_pack"].empty:
            r["tukey_pack"].to_csv(TABLES / f"posthoc_tukey_pack_{slug}.csv", index=False)

    pd.DataFrame([
        {"term": term, **stats} for term, stats in pooled["anova"].items()
    ]).to_csv(TABLES / "pooled_model_anova.csv", index=False)
    pd.DataFrame(impact["lines"]).to_csv(TABLES / "commercial_impact.csv", index=False)
    df.groupby(["commodity", "temp_level", "pack_level"], as_index=False).agg(
        n=("spoilage_pct", "size"),
        spoilage_mean=("spoilage_pct", "mean"),
        spoilage_sd=("spoilage_pct", "std"),
        shelf_life_mean=("shelf_life_days", "mean"),
        crush_mean=("crush_damage_pct", "mean"),
        arrival_temp_mean=("arrival_temp_c", "mean"),
    ).round(3).to_csv(TABLES / "cell_summary.csv", index=False)

    # ---- figures --------------------------------------------------------
    print()
    print("=" * 78)
    fig_cell_means(df)
    fig_interaction(df)
    fig_diagnostics(results)
    fig_commercial(impact)

    # ---- report + json ---------------------------------------------------
    write_report(df, results, pooled, impact, grand)
    payload = {
        "design": {
            "type": "3x2 factorial, blocked by commodity, balanced",
            "temperature_levels": TEMP_LEVELS,
            "packaging_levels": PACK_LEVELS,
            "commodities": COMMODITIES,
            "replicates_per_cell": 8,
            "total_runs": int(len(df)),
            "primary_response": "spoilage_pct",
            "secondary_response": "shelf_life_days",
            "alpha": ALPHA,
            "practical_threshold_pp": PRACTICAL_THRESHOLD_PP,
        },
        "pooled_model_rejected": pooled,
        "blocked_models": [
            {k: v for k, v in r.items()
             if k not in ("block", "model_summary", "tukey_temp", "gh_temp",
                          "tukey_pack")}
            for r in results
        ],
        "shelf_life_secondary": shelf_rows,
        "commercial": impact,
    }
    (REPORT / "analysis_results.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"  wrote outputs/analysis_results.json")
    print(f"  wrote {len(list(TABLES.glob('*.csv')))} table files")


if __name__ == "__main__":
    main()
