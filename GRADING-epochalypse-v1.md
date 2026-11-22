# Grading: calibration v1 on the Epochalypse simulated DR4 epochs

Graded on 2026-11-22, ten days before Gaia DR4. Rules: [VALIDATION-epochalypse.md](VALIDATION-epochalypse.md)
(fixed before the run). Job `rvdr4-epochalypse-v1` (repository commit c6d79bb, rvdr4.py v0.2.0), run on
2026-11-22. The 41 unit tests passed first. Runtime was 369 s and peak memory 0.38 GB. Work by AI agents
working in the name of MV2 (Munim, Inc.). No human expert has checked it yet.

Files: [results/epochalypse-v1/results_epochalypse.json](results/epochalypse-v1/results_epochalypse.json),
[systems_1_companion.csv](results/epochalypse-v1/systems_1_companion.csv) (300 systems x oracle/perturbed),
[nulls_0_companion.csv](results/epochalypse-v1/nulls_0_companion.csv) (300 nulls).

## 1. Pre-registered outcome: EPOCHALYPSE_PASS

| Gate | Sample | Value | Threshold | Result |
|---|---|---|---|---|
| G1 null false-alarm rate | 300 zero-companion systems | 0 of 300 (0.0%) | <= 1% | pass |
| G2 inclination 95% coverage (oracle, P <= 11 yr) | 221 measured systems | 200 of 221 (90.5%) | >= 85% | pass |
| G3 ORBIT_INCONSISTENT_WITH_RV rate | 226 detected systems | 5 of 226 (2.2%) | <= 10% | pass |

Reported without gates:
- **Perturbed ephemeris** (P x (1 + 0.01 N), e + 0.02 N, Tp + 0.02 P N): inclination coverage 87.3%,
  mass coverage 91.4%, same verdict mix within one system.
- **All periods up to 20 yr:** inclination coverage 86.1% (oracle) and 83.3% (perturbed); 13 of 300 not
  detected.
- **True-mass 95% coverage** (oracle, P <= 11 yr): 93.7%.
- **Median absolute inclination error:** 3.4 deg.

As pre-registered, the DR4 verdicts will be reported with this calibration attached. No correction to the
error bars is required by the v1 rules.

## 2. Diagnostic found after the run (not pre-registered)

Looking at which systems were missed, the misses are concentrated at **nearly face-on orbits**. The split at
25 deg below was chosen after seeing the data.

| True folded inclination | Measured systems | 95% interval covers truth | Median too high |
|---|---|---|---|
| < 25 deg | 21 | 14 (66.7%) | 20 of 21 |
| >= 25 deg | 200 | 186 (93.0%) | - |

For the 21 nearly face-on systems, the median of (estimated - true) inclination is **+11.1 deg**. True
masses are underestimated by the same factor, and mass coverage is also 14 of 21. Among the 23 oracle
systems (P <= 11 yr) labelled TRUE_MASS_MEASURED (that is, not confidently above 13 Jupiter masses), 12
are brown dwarfs in truth, and 6 of those 12 have true inclinations below 25 deg. The simulated population
is dominated by brown dwarfs, so this ratio says nothing about real planets. It does show the direction of
the error: a face-on brown dwarf can pass for a planet.

**Cause.** rvdr4.py reads the inclination from the shape of the fitted Thiele-Innes ellipse alone. For a
face-on orbit, noise always makes the projected ellipse look flatter, never rounder, so the inclination is
biased toward edge-on. The RV amplitude a0 sin i, which the method already computes, is used only as a
consistency check. For nearly face-on orbits that amplitude carries most of the information
(sin i = a0 sin i / a0).

**Why it matters for the pre-registered targets.** Every published DR3-based inclination for these
planets is low: GJ 463 b 21.8 deg (Morales et al. 2025) versus 27-28 deg (Sozzetti 2023; Piccinini et al.
2026), GJ 9773 b 22.3 deg, GJ 317 b/c 33.9 deg and GJ 3512 b/c 36 deg. A shape-only estimate biased toward
edge-on would favour the larger GJ 463 b inclination for the wrong reason.

**Response.** [PREREGISTRATION-amendment2.md](PREREGISTRATION-amendment2.md), dated 2026-11-22 (before
DR4), adds an RV-amplitude-constrained estimator (gaiaorb/rvincl.py). It is calibrated on the same
Epochalypse systems under pre-stated gates, and a fixed rule decides which estimator is primary on DR4 day.
The v1 verdict rules are unchanged.

## 3. Limits

- One shard of the simulator, high-SNR one-companion systems only; white Gaussian per-transit noise.
- The face-on diagnostic rests on 21 systems and was found after looking.
- The simulator's population is mostly brown dwarfs; real RV planets are lighter and their wobbles
  smaller.
