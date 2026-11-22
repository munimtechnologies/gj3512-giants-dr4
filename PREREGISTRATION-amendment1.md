# Pre-registration, Amendment 1: GJ 317 b and c, GJ 463 b and GJ 9773 b in Gaia DR4, and a test of the DR3-based inclination predictions

Frozen on 2026-11-22, ten days before Gaia DR4 (2026-12-02, 12:00 CET), before any DR4 data
exist. Version 0.2.0 of gj3512-giants-dr4. Work by AI agents working in the name of MV2
(Munim, Inc.). Apache-2.0. The original pre-registration (PREREGISTRATION.md, v0.1.0, GJ 3512)
is unchanged: same code path, same verdicts. This amendment adds targets and reported
quantities only.

## 1. Why these planets

Morales et al. 2025 (A&A; arXiv:2507.15516, "The CARMENES search for exoplanets around M dwarfs.
Revisiting the GJ 317, GJ 463, and GJ 3512 systems...") derived RV orbits for giant planets around
five M dwarfs. They estimated their inclinations and true masses from the Gaia DR3 astrometric
excess noise, adding the Hipparcos-Gaia proper-motion anomaly for GJ 463. They wrote that "these caveats will all be overcome with the
release of epoch astrometry as part of the Gaia DR4". For the two-planet systems they assumed the
orbits are coplanar.

For GJ 463 b, three published DR3-based solutions disagree:

| Source | folded inclination (deg) | true mass (M_Jup) |
|---|---|---|
| Morales et al. 2025, Table 13 | 21.8 (+1.5/-1.1) | 4.44 (+0.31/-0.29) |
| Sozzetti 2023 (A&A 670, L17), prograde / retrograde | 27 +- 3 / 28 (+3/-2) | 3.4 +- 0.3 / 3.6 +- 0.4 |
| Piccinini et al. 2026 (arXiv:2601.09401), prograde / retrograde | 27.7 (+4.6/-3.9) / 27.7 (+4.6/-3.9) | 3.3 (+0.4/-0.3) / 3.4 +- 0.4 |

Sozzetti 2023 wrote that Gaia DR4 should deliver a complete orbital solution.

**Questions.**
1. What are the inclinations and true masses of GJ 317 b, GJ 317 c, GJ 463 b and GJ 9773 b,
   measured directly from DR4 along-scan epoch astrometry with the RV elements fixed?
2. Are GJ 317 b and c (and GJ 3512 b and c) coplanar, as the DR3 analyses assumed?
3. Do the published DR3-based inclinations (excess noise, proper-motion anomaly) survive a direct
   orbit measurement? For GJ 463 b, which of the disagreeing solutions does DR4 support?

## 2. Frozen inputs (all in gaiaorb/rvdr4.py)

Stars (Gaia DR3, fetched 2026-11-22 from the Gaia archive; the DR3 excess noise of each matches
Morales et al. 2025 Table 1, confirming the identification; masses from that table):

| Star | DR3 source_id | G | parallax (mas) | RUWE | excess noise (mas) | DR3 RV (km/s) | M (Msun) |
|---|---|---|---|---|---|---|---|
| GJ 317 | 5701750715316821248 | 10.77 | 65.877 | 1.111 | 0.358 | 87.11 | 0.402 +- 0.014 |
| GJ 463 | 1584088431025057792 | 10.55 | 54.447 | 1.407 | 0.183 | 8.92 | 0.486 +- 0.014 |
| GJ 9773 | 2723557788822526720 | 10.91 | 62.878 | 1.195 | 0.133 | 8.19 | 0.358 +- 0.011 |

None has a DR3 non-single-star solution (non_single_star = 0).

Planets (Morales et al. 2025, Table 4; T_conj = time of conjunction, BJD):

| | P (d) | K (m/s) | e | omega (deg) | T_conj | m sin i (M_Jup) | a0_min (mas) | a0 at the Morales i (mas) |
|---|---|---|---|---|---|---|---|---|
| GJ 317 b | 695.69 | 71.39 | 0.0826 | 11.3 | 2458089.4 | 1.688 | **0.300** | 0.54 |
| GJ 317 c | 6500 | 29.33 | 0.217 | 106.7 | 2462735 | 1.430 | **1.127** | 2.02 |
| GJ 463 b | 3459 | 36.27 | 0.146 | 329.3 | 2454313 | 1.649 | **0.621** | 1.67 |
| GJ 9773 b | 2965 | 15.11 | 0.218 | 99.0 | 2459362 | 0.543 | **0.253** | 0.67 |

**GJ 9773 b omega.** Table 4 prints e cos(omega) = +0.034 and e sin(omega) = 0.214, which give
omega = 81.0 deg, but it also prints omega = 99.0 deg. That sign conflict moves the periastron
time by about 90 days. The primary fit uses 99.0 deg; the identical fit with 81.0 deg is also run.

GJ 508.2 b, also in that paper, has only lower limits on its period (> 5300 d), so it is not
tested.

## 3. Data

As in PREREGISTRATION.md, Section 3. Gaia DR4 epoch astrometry from the archive DataLink is
collapsed to per-transit along-scan normal points by gaiaorb/epochs.py. Only the I/O may be adapted
after release.

## 4. Analysis (frozen; gaiaorb/rvdr4.py v0.2.0, rv_job_a1.py)

Identical to PREREGISTRATION.md, Section 4, per system:
- GJ 317: b and c fitted jointly (5 + 8 parameters).
- GJ 463 and GJ 9773: one planet each (5 + 4).

The DR3 RV fixes the perspective acceleration. Additions, all report-only unless stated otherwise:
1. The posterior of the folded inclination, min(i, 180 - i), from the same 4,000 draws.
2. Mutual inclination of two-planet systems (GJ 317, GJ 3512): 4,000 joint draws of the eight TI
   constants (seed 20261203). Astrometry alone cannot distinguish Omega from Omega + 180 deg, so
   both branches are reported ("same", "flipped"), plus the branch picked by matching each planet's
   astrometric omega to its RV omega (convention-dependent, flagged as such).
3. GJ 9773 b: the omega = 81.0 deg fit, reported in full (chi^2, jitter, elements, posterior).

## 5. Verdicts

- **Per planet:** the verdicts of PREREGISTRATION.md, Section 5 (TRUE_MASS_MEASURED,
  TRUE_MASS_MEASURED_BROWN_DWARF, ORBIT_INCONSISTENT_WITH_RV, RV_ORBIT_NOT_SEEN, NOT_DETECTED,
  INCONCLUSIVE). They require the shared DR4 validation (proxima-dr4-validate) to pass.
- **GJ 9773 b:** the verdict stands only if the omega = 99.0 and 81.0 deg fits give the same
  verdict; otherwise **INCONCLUSIVE_OMEGA_AMBIGUOUS**.
- **Prediction checks (secondary, pre-registered):** applied for each planet with a
  TRUE_MASS_MEASURED(_BROWN_DWARF) verdict and each published prediction in rvdr4.PREDICTIONS.
  These are Morales et al. 2025 Table 13 for GJ 3512 b/c, GJ 317 b/c, GJ 463 b and GJ 9773 b,
  plus Sozzetti 2023 and Piccinini et al. 2026 for GJ 463 b.
  - **AGREES:** the predicted folded inclination lies inside the DR4 95% interval (2.5th-97.5th
    percentile) of the folded inclination.
  - **DISAGREES:** otherwise.
  - **NOT_EVALUATED:** the planet's verdict is not a measurement.
  - Also reported: the offset in units of the DR4 68% half-width, and whether the predicted true
    mass lies inside the DR4 95% true-mass interval.
- **GJ 463 b reading, fixed now:**
  - If Morales et al. AGREES and both other papers DISAGREE (or the reverse), DR4 favours that
    solution.
  - If all agree, DR4 cannot separate them.
  - If none agrees, all three DR3-based solutions failed.

## 6. What we have already seen (disclosure)

The papers above, and the Gaia DR3 catalogue values for the three stars. No DR4 data exist. The
code has been tested only on synthetic epochs (test_rvdr4_a1.py).

## 7. Limits stated in advance

- GJ 317 c (P = 6500 d) has a period about three times DR4's 66-month baseline. Much of its signal
  is absorbed by the proper motion, and its inclination may be weakly constrained; the covariance
  carries this. The RV period and phase errors of c (+-169 d, +-179 d) are not propagated into the
  primary fit.
- GJ 463 b (P = 3459 d) and GJ 9773 b (P = 2965 d) also exceed the baseline, but less so.
- These stars are bright (G = 10.6-10.9) and their signals large. DPAC's own DR4 pipeline may
  publish orbits for some of them on release day. Ours is an independent, pre-registered
  measurement with the RV elements fixed and multi-planet systems fitted jointly. Any DPAC
  solution will be compared, not used.
- Undetected further companions would bias the fits.
