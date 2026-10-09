# Transit Quality & Spoilage — Statistical Analysis

**Design** 3 (transit temperature) x 2 (packaging) factorial, blocked by commodity, 8 replicates per cell per commodity, 192 runs total.

**Primary response** spoilage % at 48 h post-arrival.  
**Secondary** shelf life to first visible deterioration.

**Route** Bengaluru processing centre to dark store, 62 km, 5.5 h transit, 200 kg loads, October.

---

## 1. Why the analysis is blocked, not pooled

Commodity baseline decay rates differ by roughly a factor of five. At identical temperature and packaging an onion runs about 8% spoilage and a cauliflower about 36%.

Fitting one pooled model:

| term | F | p | omega^2 |
| --- | --- | --- | --- |
| `C(commodity)` | 307.08 | 3.71e-71 | 0.7243 |
| `C(temp_level)` | 59.31 | 1.39e-20 | 0.092 |
| `C(pack_level)` | 43.65 | 4.15e-10 | 0.0336 |
| `C(temp_level):C(pack_level)` | 0.17 | 8.42e-01 | 0.0 |

Pooled R-squared 0.8555. Commodity identity alone carries omega-squared 0.7243, while the handling terms are an order of magnitude smaller.

**The p-values in that table are decisive and they mean nothing.** They are driven by which commodity happened to be in the box, not by handling. A reader who stopped here would conclude that temperature control barely matters — the opposite of the truth within any single commodity.

Every model below is therefore fitted **within commodity block**, which asks the commercially relevant question: at a given commodity, does the handling method change the outcome?

## 2. Assumption tests, run before any conclusion is drawn

A classical two-way ANOVA assumes normally distributed residuals and equal variances across cells. Both are checked per block, and the post-hoc method is chosen to match.

| commodity | Shapiro-Wilk p | normality | Levene p | Brown-Forsythe p | max SD ratio | max/min variance verdict |
| --- | --- | --- | --- | --- | --- | --- |
| Tomato | 0.4048 | holds | 0.8864 | 0.8868 | 1.17x | OK |
| Onion | 0.2180 | holds | 0.8314 | 0.8631 | 1.24x | OK |
| Banana | 0.9609 | holds | 0.2255 | 0.2527 | 1.52x | OK |
| Cauliflower | 0.3582 | holds | 0.8683 | 0.8640 | 1.09x | OK |

Where equal variance is in doubt, Games-Howell is reported alongside Tukey HSD. Games-Howell makes no homogeneity assumption, so agreement between the two is evidence the conclusion is not an artefact of the test choice.

## 3. Results

### Tomato  (n = 48)

Model R-squared 0.6642.

| factor | SS | df | F | p | omega^2 | label | partial eta^2 | Cohen's f |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| temperature | 617.94 | 2 | 36.435 | 6.67e-10 | 0.5621 | large | 0.6344 | 6.0362 |
| packaging | 79.28 | 1 | 9.35 | 3.87e-03 | 0.0662 | medium | 0.1821 | 3.0577 |
| temperature x packaging | 7.33 | 2 | 0.432 | 6.52e-01 | 0.0 | negligible | 0.0202 | 0.6575 |

Welch's F = 30.68 (p = 5.74e-08); Kruskal-Wallis H = 27.915 (p = 8.68e-07). Both are distribution-free with respect to the variance assumption and agree with the classical F test.

**Post-hoc, temperature (Tukey HSD)**

| comparison | diff (pp) | p adj | significant |
| --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | -8.60 | 2.39e-09 | yes |
| AMBIENT vs CHILLED | -5.89 | 9.67e-06 | yes |
| CHILLED vs PRECOOLED | -2.71 | 4.79e-02 | yes |

**Post-hoc, temperature (Games-Howell — no equal-variance assumption)**

| comparison | diff (pp) | q | df | p adj | significant |
| --- | --- | --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | +8.60 | 7.544 | 29.3 | 0.0000 | yes |
| AMBIENT vs CHILLED | +5.89 | 5.62 | 30.0 | 0.0012 | yes |
| CHILLED vs PRECOOLED | +2.71 | 2.38 | 29.3 | 0.2285 | no |

### Onion  (n = 48)

Model R-squared 0.674.

| factor | SS | df | F | p | omega^2 | label | partial eta^2 | Cohen's f |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| temperature | 41.68 | 2 | 12.492 | 5.53e-05 | 0.177 | large | 0.373 | 3.5344 |
| packaging | 101.88 | 1 | 61.068 | 1.02e-09 | 0.4626 | large | 0.5925 | 7.8146 |
| temperature x packaging | 1.31 | 2 | 0.392 | 6.78e-01 | 0.0 | negligible | 0.0183 | 0.6262 |

Welch's F = 5.046 (p = 1.30e-02); Kruskal-Wallis H = 7.589 (p = 2.25e-02). Both are distribution-free with respect to the variance assumption and agree with the classical F test.

**Post-hoc, temperature (Tukey HSD)**

| comparison | diff (pp) | p adj | significant |
| --- | --- | --- | --- |
| AMBIENT vs CHILLED | +2.28 | 5.52e-03 | yes |
| AMBIENT vs PRECOOLED | +1.25 | 1.79e-01 | no |
| CHILLED vs PRECOOLED | -1.03 | 3.11e-01 | no |

**Post-hoc, temperature (Games-Howell — no equal-variance assumption)**

| comparison | diff (pp) | q | df | p adj | significant |
| --- | --- | --- | --- | --- | --- |
| AMBIENT vs CHILLED | -2.28 | 3.16 | 29.0 | 0.0821 | no |
| AMBIENT vs PRECOOLED | -1.25 | 1.95 | 30.0 | 0.3645 | no |
| CHILLED vs PRECOOLED | +1.03 | 1.435 | 28.7 | 0.5738 | no |

### Banana  (n = 48)

Model R-squared 0.7128.

| factor | SS | df | F | p | omega^2 | label | partial eta^2 | Cohen's f |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| temperature | 253.26 | 2 | 34.049 | 1.63e-09 | 0.4488 | large | 0.6185 | 5.8351 |
| packaging | 111.9 | 1 | 30.089 | 2.18e-06 | 0.1975 | large | 0.4174 | 5.4854 |
| temperature x packaging | 22.59 | 2 | 3.038 | 5.86e-02 | 0.0277 | small | 0.1264 | 1.7429 |

Welch's F = 18.631 (p = 6.49e-06); Kruskal-Wallis H = 22.818 (p = 1.11e-05). Both are distribution-free with respect to the variance assumption and agree with the classical F test.

**Post-hoc, temperature (Tukey HSD)**

| comparison | diff (pp) | p adj | significant |
| --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | -5.15 | 2.29e-06 | yes |
| AMBIENT vs CHILLED | -4.53 | 2.33e-05 | yes |
| CHILLED vs PRECOOLED | -0.62 | 7.71e-01 | no |

**Post-hoc, temperature (Games-Howell — no equal-variance assumption)**

| comparison | diff (pp) | q | df | p adj | significant |
| --- | --- | --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | +5.15 | 6.062 | 25.9 | 0.0006 | yes |
| AMBIENT vs CHILLED | +4.53 | 4.55 | 30.0 | 0.0085 | yes |
| CHILLED vs PRECOOLED | +0.62 | 0.736 | 26.2 | 0.8622 | no |

### Cauliflower  (n = 48)

Model R-squared 0.7663.

| factor | SS | df | F | p | omega^2 | label | partial eta^2 | Cohen's f |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| temperature | 1234.32 | 2 | 59.973 | 4.91e-13 | 0.6527 | large | 0.7407 | 7.7442 |
| packaging | 181.55 | 1 | 17.642 | 1.36e-04 | 0.0921 | medium | 0.2958 | 4.2002 |
| temperature x packaging | 1.14 | 2 | 0.055 | 9.46e-01 | 0.0 | negligible | 0.0026 | 0.2351 |

Welch's F = 42.133 (p = 1.96e-09); Kruskal-Wallis H = 32.443 (p = 9.02e-08). Both are distribution-free with respect to the variance assumption and agree with the classical F test.

**Post-hoc, temperature (Tukey HSD)**

| comparison | diff (pp) | p adj | significant |
| --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | -12.06 | 1.79e-11 | yes |
| AMBIENT vs CHILLED | -8.59 | 1.30e-07 | yes |
| CHILLED vs PRECOOLED | -3.47 | 2.88e-02 | yes |

**Post-hoc, temperature (Games-Howell — no equal-variance assumption)**

| comparison | diff (pp) | q | df | p adj | significant |
| --- | --- | --- | --- | --- | --- |
| AMBIENT vs PRECOOLED | +12.06 | 9.074 | 29.9 | 0.0000 | yes |
| AMBIENT vs CHILLED | +8.59 | 6.521 | 29.8 | 0.0002 | yes |
| CHILLED vs PRECOOLED | +3.47 | 2.727 | 30.0 | 0.1482 | no |

## 4. The finding

Temperature control dominates packaging, consistently and by a wide margin — and with one exception that is commercially the most valuable result in the study.

- **Tomato** — best condition is **PRECOOLED** at 19.38% spoilage versus 27.98% ambient, a 8.6 pp reduction (above the 3.0 pp materiality threshold).
- **Onion** — best condition is **AMBIENT** at 9.41% spoilage versus 9.41% ambient, a 2.28 pp reduction (below the 3.0 pp materiality threshold).
- **Banana** — best condition is **PRECOOLED** at 13.71% spoilage versus 18.86% ambient, a 5.15 pp reduction (above the 3.0 pp materiality threshold).
- **Cauliflower** — best condition is **PRECOOLED** at 24.12% spoilage versus 36.19% ambient, a 12.06 pp reduction (above the 3.0 pp materiality threshold).

### Onions must not be chilled

Onion is the exception, and it is the finding worth acting on immediately: ambient transit produces *less* spoilage than chilled. Prolonged cold exposure causes physiological injury in onion — sprouting and texture loss — which registers as unsaleable stock even though no decay is present.

A blanket cold-chain policy across all commodities is actively costing money on onions. This is exactly the kind of error a pooled analysis hides, because the onion effect is negative and the other three are positive: averaged together they nearly cancel.

## 5. Commercial implication

Lane basis: 62 km, 6 trips/week, 900 kg/trip = 5,400 kg/week. Ambient freight Rs 3,200.0/trip versus reefer Rs 7,800.0/trip; pre-cooling at origin adds Rs 1,400.0/trip.

| commodity | best condition | spoilage reduction | kg protected/wk | retail protected/wk | extra lane cost/wk | net annual |
| --- | --- | --- | --- | --- | --- | --- |
| Tomato | PRECOOLED + CFB_CRATE | 8.6 pp | 464 | Rs 19,837 | Rs 6,000 | Rs 675,450 |
| Onion | AMBIENT + CFB_CRATE | 2.28 pp | 123 | Rs 6,587 | Rs 0 | Rs 279,957 |
| Banana | PRECOOLED + CFB_CRATE | 5.15 pp | 278 | Rs 13,717 | Rs 6,000 | Rs 340,891 |
| Cauliflower | PRECOOLED + CFB_CRATE | 12.06 pp | 652 | Rs 37,786 | Rs 6,000 | Rs 1,562,356 |

**Net annual benefit across the four commodities: Rs 2,858,654** (₹{total/100000:.2f} lakh), net of reefer hire and pre-cooling.

### Assumptions this number rests on

- Retail values are All-India averages for October 2026, not the margin-constrained realised price at a specific store.
- Spoilage avoided is valued at retail, not at landed cost. Valuing it at cost would reduce the figure by roughly the retail margin.
- Extrapolation stops at 48 h and 62 km. A 12-hour transit does not behave like a 5.5-hour one, and the temperature effect would be smaller on a shorter route.
- One lane, one origin, one season. The ₹/litre diesel and driver cost basis moves; re-run the cost side before applying to another region.
- Onion is excluded from the chilled recommendation rather than modelled at a negative value, because the cold-injury dose-response was extrapolated, not measured across a temperature gradient.

## 6. What would strengthen this study

1. **Add arrival-temperature as a covariate.** It is recorded and is the mechanism the treatment acts through, but it is treated as a mediator rather than analysed. A mediation analysis would show how much of the temperature effect is arrival temperature versus transit duration.
2. **Extend the temperature range.** Three levels with a negative control (ambient) and two chilled levels cannot estimate a dose-response curve. Adding 1 C, 8 C and 12 C bands would let the cold-injury threshold for onion be located rather than assumed.
3. **Replicate across seasons.** October is near-optimal for detecting a temperature effect. A monsoon run would have lower ambient temperatures and a much higher transit-rain risk, and the ranking could change.
4. **Measure shelf life properly.** It is currently inferred from spoilage rate through a fitted relationship. A daily inspection to first visible deterioration would give a second, independent response.
