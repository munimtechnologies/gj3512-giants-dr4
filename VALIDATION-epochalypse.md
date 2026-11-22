# Calibration of the DR4 fitter on simulated Gaia DR4 epoch astrometry (Epochalypse)

Fixed on 2026-11-22, before the job runs and ten days before Gaia DR4. Version 0.1.0. Work by AI
agents working in the name of MV2 (Munim, Inc.). Apache-2.0.

## 1. Why

Every Gaia DR4 verdict in this repository (GJ 3512, GJ 317, GJ 463, GJ 9773) comes from one frozen
method, gaiaorb/rvdr4.py. It is a linear Thiele-Innes fit of the along-scan epochs, with the period,
eccentricity and periastron time fixed from radial velocities, and a jitter term. Until now it has
been tested only on its own synthetic data. Yahalomi et al. 2026 (arXiv:2611.06915; Zenodo
10.5281/zenodo.23789041, CC BY 4.0) released the "Epochalypse" catalogue:
- DR4-like per-transit epoch astrometry for millions of real Gaia stars, with injected companions
  of known orbits.
- The scanning law follows GOST.
- The noise is UEVA-based, with a smaller reported formal error, as in real Gaia data.
- Matched zero-companion systems serve as nulls.

The paper benchmarks only a periodogram; it reports no orbit-fitting pipelines. Here the frozen
method is run, unchanged, on this independent simulator, so its error bars and false-alarm rate can
be checked before DR4.

## 2. Convention independence

The four Thiele-Innes constants span every linear map from the orbit's (X, Y) to the sky offsets.
The fitted ellipse, its chi^2 improvement, the reflex semi-major axis a0 and the folded inclination
min(i, 180 - i) therefore do not depend on the simulator's sign conventions for the scan angle,
node or argument of periastron. Only these convention-free quantities are graded.

## 3. Data and selection (frozen in epocha_job.py)

- Truth: `injected_solutions_1_companion_high_snr.parquet` (all one-companion systems with
  SNR_tot >= 5).
- Epochs: the first `epochs_rank*.parquet` shard streamed from `simulated_astrometry_1_companion.tar`.
  Nulls: the first shard of `simulated_astrometry_0_companion.tar`.
- Systems: one-companion systems in that shard with:
  - 0.5 <= P <= 20 yr (P from Kepler's law with the truth masses);
  - reflex semi-major axis 0.1 <= alpha <= 3 mas (alpha = M_p/(M_star+M_p) a parallax);
  - e <= 0.6;
  - at least 40 transits.
  Sorted by source id; the first 300.
- Nulls: the first 300 zero-companion systems (sorted by source id) with at least 40 transits.
  Each is given the ephemeris of one selected one-companion system, in turn.

## 4. What is run

For each system, `rvdr4.run` and `rvdr4.verdict` (validated = True), exactly as on DR4. The star has
no perspective term (the simulator has none). Inputs:
- Period, eccentricity and periastron time from the truth.
- K chosen so that a0_min = alpha sin i.
- m sin i = M_p sin i (+-2%), stellar mass +-2%.

Two modes:
- **oracle:** the truth ephemeris.
- **perturbed:** P x (1 + 0.01 N), e + 0.02 N, T_p + 0.02 P N (seed 20261122), mimicking RV
  ephemeris errors.

## 5. Gates (oracle mode, P <= 11 yr = two DR4 baselines)

- **G1 null:** fraction of null fits with dchi2 >= 23.51 is <= 1% (the nominal rate is 1e-4).
  Needs >= 100 nulls.
- **G2 inclination coverage:** among TRUE_MASS_MEASURED(_BROWN_DWARF) systems, the true folded
  inclination lies inside the 95% interval in >= 85% of cases. Needs >= 30 systems.
- **G3 consistency:** among systems with dchi2 >= 23.51, the fraction with
  ORBIT_INCONSISTENT_WITH_RV is <= 10%. Needs >= 30 systems.

**Outcome:** EPOCHALYPSE_PASS if all three pass, EPOCHALYPSE_FAIL:<gates> otherwise, and
EPOCHALYPSE_UNDERPOWERED if any sample is too small.

Reported without a gate:
- the same numbers in perturbed mode and for all periods up to 20 yr;
- true-mass coverage;
- the verdict mix;
- the median inclination error.

## 6. What follows

- **PASS:** the DR4 verdicts are reported with this calibration attached.
- **FAIL:** before DR4 opens (2026-12-02), a dated amendment states how the DR4 error bars will be
  corrected, for example inflating the TI covariance by the measured coverage deficit, or that the
  affected verdicts will be reported as provisional. The FAIL itself is reported.

## 7. Limits

- The simulator's noise is white and Gaussian per transit. Real DR4 data have CCD-level structure,
  attitude noise and outliers.
- Only one shard of each population is used. These are the high-SNR one-companion systems; no
  two-companion systems are used, so the coplanarity statistic is not calibrated here.
- The perturbed mode's ephemeris errors are generic, not those of the real RV fits.
