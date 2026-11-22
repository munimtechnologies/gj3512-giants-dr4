# Pre-registration: the true masses of GJ 3512 b and c from Gaia DR4 epoch astrometry

Frozen on 2026-11-22, ten days before Gaia DR4 (2026-12-02, 12:00 CET), before any DR4 data
exist. Version 0.1.0. Work by AI agents working in the name of MV2 (Munim, Inc.). Apache-2.0.

## 1. The question

GJ 3512 is an M5.5 dwarf of 0.123 solar masses at 9.5 pc. Morales et al. 2019 (Science 365,
1441; arXiv:1909.12174) found a giant planet (m sin i = 0.46 M_Jup, P = 204 d, e = 0.43) around it,
a planet that core accretion struggles to form around so small a star; they suggested disc
instability. Morales et al. 2025 (A&A; arXiv:2507.15516) confirmed a second giant, c (m sin i =
0.45 M_Jup, P = 2354 d). Only minimum masses are known. Morales et al. 2025 limited the true
masses with the Gaia DR3 excess noise (0.264 mas) alone: about 0.78 (+0.04/-0.27) M_Jup for both,
at i about 36 (+27/-1.2) deg, Omega about 142.5 deg, assuming coplanar orbits. They wrote that "these caveats will all be
overcome with the release of epoch astrometry as part of the Gaia DR4". Morales et al. 2019 had
already forecast a reflex motion of about 120 microarcseconds, twice Gaia's accuracy.

**Question.** What are the inclinations and true masses of GJ 3512 b and c, measured directly
from the Gaia DR4 along-scan epoch astrometry with the RV orbital elements fixed? Are the two
orbits coplanar? Is either companion heavier than 13 M_Jup?

## 2. Frozen inputs

- Star: Gaia DR3 1040827088726734976 (G 13.08, parallax 105.294 mas, RUWE 1.25,
  non_single_star 0, RV 7.13 km/s), mass 0.123 +- 0.009 Msun (Morales et al. 2025).
- Planets: Morales et al. 2025, Table 4 (all in gaiaorb/rvdr4.py):

| | P (d) | K (m/s) | e | omega (deg) | T_conj (BJD) | m sin i (M_Jup) | a0_min (mas) |
|---|---|---|---|---|---|---|---|
| b | 203.109 | 71.28 | 0.4279 | 126.80 | 2457739.51 | 0.461 | **0.1266** |
| c | 2354 | 28.02 | 0.089 | 50 | 2457339 | 0.453 | **0.6359** |

  a0_min = K sqrt(1-e^2) P / (2 pi) x parallax is the star's reflex semi-major axis for
  sin i = 1; Gaia must see a0 = a0_min / sin i. Periastron times follow from T_conj with the
  true anomaly 90 deg - omega at conjunction (radvel/juliet convention).

## 3. Data

Gaia DR4 epoch astrometry for the star, fetched from the Gaia archive DataLink and collapsed to
one along-scan normal point per field-of-view transit by the shared gaiaorb code
(gaiaorb/epochs.py, as in proxima-c-dr4, radio-planets-dr4 and pceb-planets-dr4). Only the I/O
may be adapted after the release (source lookup, URLs, file format); nothing in Sections 4-5.

## 4. Analysis (frozen; gaiaorb/rvdr4.py)

1. Model: 5 single-star parameters (exact perspective acceleration with the DR3 RV) + 4
   Thiele-Innes constants per planet with P, e and periastron time fixed at the RV values
   (13 linear parameters, one exact weighted least-squares solve). One jitter term makes the
   reduced chi^2 of this full model 1; one pass of 5-sigma clipping against the full model.
2. For each planet: dchi2 = chi^2(model without that planet) - chi^2(full), same errors.
3. Campbell elements (a0, i, omega, Omega) from the TI constants; posteriors by 4,000 Monte
   Carlo draws from the TI covariance (seed 20261202), with m sin i and the stellar mass drawn
   from their errors. True mass = m sin i / sin i; cross-check: mass from a0 alone.
4. Consistency with the RVs: z = |median(a0 sin i) - a0_min| / half-width of its 68% interval.
5. Expected dchi2 of the minimum-mass orbit (i = 90 deg, median over 12 node angles) with the
   real epochs and errors.

## 5. Verdict per planet (requires the shared DR4 validation, proxima-dr4-validate, to pass)

- **TRUE_MASS_MEASURED**: dchi2 >= 23.51 (4 dof, p = 1e-4) and z <= 3.
- **TRUE_MASS_MEASURED_BROWN_DWARF**: as above and the 16th percentile of the true mass > 13 M_Jup.
- **ORBIT_INCONSISTENT_WITH_RV**: dchi2 >= 23.51 but z > 3.
- **RV_ORBIT_NOT_SEEN**: dchi2 < 23.51 although the minimum-mass orbit should give >= 50
  (this would contradict the RV planet or the DR4 errors).
- **NOT_DETECTED**: dchi2 < 23.51 and the minimum-mass orbit would give < 50.
- **INCONCLUSIVE**: validation failed, no epochs, or an error.

Reported without a verdict: inclinations of b and c and their difference (coplanarity),
omega(Gaia) vs omega(RV) modulo 180 deg, and agreement with the excess-noise masses
(0.78 +0.04/-0.27 M_Jup) and the DR3-based inclination range of Morales et al. 2025.

## 6. What we have already seen (disclosure)

The papers above; Gaia DR3 catalogue values for the star. No DR4 data exist yet. Code tested
only on synthetic epochs.

## 7. Limits stated in advance

- Fixed RV elements: errors in P, e, periastron time of c (P 2354 d is longer than DR4's
  66-month baseline) are not propagated into the primary fit; c's orbit is partly absorbed by
  proper motion, which the covariance reflects.
- The Thiele-Innes fit leaves omega free, so the RV omega serves as a check, not a constraint.
- A third, unseen companion would bias the fit (Morales et al. 2025 exclude Saturn-mass
  planets inside ~40 yr).
- DPAC's own DR4 exoplanet pipeline may also publish an orbit for this star; our fit is an
  independent, pre-registered measurement with both RV planets fitted jointly.
