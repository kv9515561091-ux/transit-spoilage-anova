"""Generate the transit handling experiment dataset.

Design
------
A 3 x 2 factorial, blocked by commodity, fully replicated.

    Factor A: transit temperature control   3 levels
      - AMBIENT      open body, ~32 C at the tail of the route
      - CHILLED      reefer 2-6 C, loaded at origin
      - PRECOOLED    6 h pre-cool at origin, then reefer 2-4 C

    Factor B: packaging                     2 levels
      - JUTE_BAG     breathable woven bag, no stack protection
      - CFB_CRATE    corrugated ventilated crate, stackable, rigid walls

    Blocking factor: commodity              4 levels
      - Tomato, Onion, Banana, Cauliflower

    Replicates per cell per commodity:      8

    Total runs: 3 x 2 x 4 x 8 = 192

Responses
--------
    spoilage_pct          primary. Share of load unsaleable at 48 h post-arrival.
    shelf_life_days       secondary. Days to first visible deterioration.
    arrival_temp_c        process variable, recorded not controlled -- it is
                          what the temperature factor actually *does* to the load.
    crush_damage_pct      covariate. Bruising is a packaging mechanism, not a
                          temperature one, so it is separated rather than
                          left in the residual.
    moisture_loss_pct     covariate. Chilled transit of unwashed produce
                          condenses; this is the mechanism by which
                          "chilled" can lose to "ambient".

Why blocked rather than pooled
------------------------------
Commodities have wildly different baseline decay rates -- an onion at 8% and
a cauliflower at 40% at the same temperature and packaging. Pooling them into
one ANOVA would produce a significant result driven by commodity identity
while saying nothing about handling. Blocking on commodity and analysing
within blocks asks the question actually of interest: at a given commodity,
does the handling method change the outcome?
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "transit_experiment.csv"

SEED = 20261016

# --------------------------------------------------------------------------
# Factor levels
# --------------------------------------------------------------------------
TEMPERATURE_LEVELS = ["AMBIENT", "CHILLED", "PRECOOLED"]
PACKAGING_LEVELS = ["JUTE_BAG", "CFB_CRATE"]
COMMODITIES = ["Tomato", "Onion", "Banana", "Cauliflower"]
REPLICATES = 8

# --------------------------------------------------------------------------
# Per-commodity baseline spoilage at 48 h, under the WORST cell
# (AMBIENT + JUTE). Commercially anchored: transit-stage loss in the Indian
# F&V chain is widely reported at 5-6%, and these are the rates a 48 h
# ambient run in October heat would actually produce.
# --------------------------------------------------------------------------
COMMODITY_BASE = {
    #        ambient+jute baseline,   residual sd,   shelf life base (days)
    "Tomato":       (28.0, 6.5, 7.0),
    "Onion":        (7.5, 2.2, 28.0),
    "Banana":       (18.0, 5.0, 8.0),
    "Cauliflower":  (36.0, 8.0, 6.0),
}

# Temperature effect: percentage reduction in spoilable mass, and the
# arrival temperature the vehicle actually delivers.
#
# PRECOOLED beats CHILLED on arrival temperature because a reefer loaded with
# field-warm produce spends its first hours pulling the load down rather than
# holding it -- the load's own heat is part of the transit load.
TEMPERATURE_EFFECT = {
    #        spoilage reduction,  arrival temp mean, arrival temp sd
    "AMBIENT":    (0.00, 30.8, 2.4),
    "CHILLED":    (0.34, 7.2, 1.6),
    "PRECOOLED":  (0.44, 4.6, 1.1),
}

# Packaging effect on spoilage is deliberately modest. Packaging's real value
# shows up in crush damage and in stackability, not in decay rate -- and
# claiming otherwise would be overstating what the data supports.
PACKAGING_EFFECT = {
    #        spoilage reduction, crush mean, crush sd
    "JUTE_BAG":  (0.00, 7.8, 2.6),
    "CFB_CRATE": (0.09, 3.1, 1.4),
}

# Condensation penalty. Chilling unwashed produce pulls water out of the air
# onto the surface, which promotes microbial growth. This is why PRECOOLED is
# not simply "CHILLED, more so" -- it carries a slightly higher wet-surface
# penalty, and it is the mechanism that stops chilled from winning outright.
CONDENSATION_PENALTY_PCT = {
    "AMBIENT": 0.0,
    "CHILLED": 1.9,
    "PRECOOLED": 2.6,
}

MOISTURE_LOSS = {
    "AMBIENT": (3.1, 1.0),
    "CHILLED": (1.4, 0.6),
    "PRECOOLED": (1.0, 0.5),
}

# Commodities that tolerate chilling. Onions are damaged by prolonged cold
# exposure, which is a real and commercially important exception to "chilled
# is better".
CHILL_SENSITIVE = {"Onion"}

ROUTE = "Bengaluru PC -> dark store, 62 km"
TRANSIT_HOURS = 5.5
LOAD_SIZE_KG = 200.0

# Fixed covariates held constant by design, recorded for auditability.
DISPATCH_HOUR = "04:30"
ORIGIN = "Bengaluru Rural APMC"


def clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def build_rows(seed: int = SEED) -> list[dict]:
    rng = random.Random(seed)
    rows: list[dict] = []
    run_no = 0

    for commodity in COMMODITIES:
        base_pct, base_sd, base_shelf = COMMODITY_BASE[commodity]
        for temp_level in TEMPERATURE_LEVELS:
            for pack_level in PACKAGING_LEVELS:
                temp_red, temp_mean, temp_sd = TEMPERATURE_EFFECT[temp_level]
                pack_red, crush_mean, crush_sd = PACKAGING_EFFECT[pack_level]
                cond_penalty = CONDENSATION_PENALTY_PCT[temp_level]
                ml_mean, ml_sd = MOISTURE_LOSS[temp_level]

                # Cold sensitivity reverses the sign of chilling for onions.
                chill_sign = -1 if commodity in CHILL_SENSITIVE else 1

                for rep in range(1, REPLICATES + 1):
                    run_no += 1

                    arrival_temp = max(1.8, rng.gauss(temp_mean, temp_sd))

                    # Spoilage: baseline scaled by both effects, then the
                    # condensation penalty and a cold-injury term applied.
                    reduction = 1 - (1 - temp_red) * (1 - pack_red)
                    expected = base_pct * (1 - reduction)

                    if chill_sign < 0 and temp_level != "AMBIENT":
                        # Prolonged chilling of onions: physiological injury
                        # showing up as sprouting and texture loss, counted
                        # as unsaleable.
                        cold_injury = (base_pct * 0.30
                                       * (arrival_temp - 2.0) / 4.0)
                        expected += cold_injury

                    expected += cond_penalty

                    # Crush damage is a packaging mechanism and is largely
                    # independent of temperature, so it is simulated
                    # separately rather than folded into the spoilage noise.
                    crush = clamp(rng.gauss(crush_mean, crush_sd), 0.2, 30.0)
                    # Bruising is itself a spoilage route: damaged tissue
                    # respires faster and rots first.
                    expected += crush * 0.42

                    moisture = clamp(rng.gauss(ml_mean, ml_sd), 0.1, 12.0)
                    # Weight loss in transit is not spoilage -- it is
                    # dehydration. Kept as its own response for that reason.
                    expected -= moisture * 0.12

                    spoilage = clamp(rng.gauss(expected, base_sd * 0.42),
                                     0.1, 92.0)

                    # Shelf life: a logarithmic relationship with spoilage
                    # rate, which is why a percentage reduction in spoilage
                    # translates into more than a proportional gain in days.
                    shelf = base_shelf * (1 - 0.55 * spoilage / 100.0)
                    shelf = clamp(rng.gauss(shelf, base_shelf * 0.11), 0.5, 60.0)

                    # Fill rate is a downstream consequence, not a response
                    # to the handling treatment, but it is what the buyer
                    # actually feels.
                    fill_rate = clamp(100 - spoilage * 0.62, 20.0, 99.5)

                    rows.append({
                        "run_id": f"R{run_no:03d}",
                        "commodity": commodity,
                        "temp_level": temp_level,
                        "pack_level": pack_level,
                        "replicate": rep,
                        "route": ROUTE,
                        "transit_hours": TRANSIT_HOURS,
                        "load_size_kg": LOAD_SIZE_KG,
                        "dispatch_hour": DISPATCH_HOUR,
                        "origin": ORIGIN,
                        "arrival_temp_c": round(arrival_temp, 1),
                        "spoilage_pct": round(spoilage, 2),
                        "shelf_life_days": round(shelf, 2),
                        "crush_damage_pct": round(crush, 2),
                        "moisture_loss_pct": round(moisture, 2),
                        "fill_rate_pct": round(fill_rate, 2),
                        "unsaleable_kg": round(LOAD_SIZE_KG * spoilage / 100.0, 2),
                    })

    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate transit handling experiment data.")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    rows = build_rows(args.seed)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    headers = list(rows[0].keys())
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=headers)
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    print(f"[generate_experiment] seed={args.seed}  seed={SEED}")
    print(f"  wrote {OUT.relative_to(ROOT)}  ({n} runs)")
    print(f"  design: 3 temperature x 2 packaging x {len(COMMODITIES)} commodities "
          f"x {REPLICATES} replicates = {3*2*len(COMMODITIES)*REPLICATES}")
    print(f"  expected {n} rows -> {'OK' if n == 3*2*len(COMMODITIES)*REPLICATES else 'MISMATCH'}")

    # Cell means, so the data can be sanity-checked against the design.
    print()
    print("  cell means (spoilage %):")
    for commodity in COMMODITIES:
        line = [f"    {commodity:14s}"]
        for t in TEMPERATURE_LEVELS:
            for p in PACKAGING_LEVELS:
                vals = [r["spoilage_pct"] for r in rows
                        if r["commodity"] == commodity and r["temp_level"] == t
                        and r["pack_level"] == p]
                line.append(f"{t[:4]}/{p[:4]}={sum(vals)/len(vals):5.1f}")
        print("  ".join(line))


if __name__ == "__main__":
    main()
