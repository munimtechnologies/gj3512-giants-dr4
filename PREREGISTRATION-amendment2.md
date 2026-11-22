# Amendment 2: an RV-amplitude-constrained inclination, calibrated before Gaia DR4

Fixed on **2026-11-22**. That is before Gaia DR4 (2026-12-02) and before the calibration job
`rvdr4-epochalypse-ce-v1` has run. Work by AI agents working in the name of MV2 (Munim, Inc.). Apache-2.0
(code), CC BY 4.0 (text). No human expert has checked it yet.

## 1. Why

The pre-registered v1 calibration passed (see [GRADING-epochalypse-v1.md](GRADING-epochalypse-v1.md)).
A diagnostic done afterwards, and **not pre-registered**, found that the shape-only inclination of
rvdr4.py v0.2.0 is biased toward edge-on for nearly face-on orbits. Below 25 deg true inclination, 14 of 21
intervals covered the truth, and the median was too high in 20 of 21 (median error +11 deg). All six
target planets have published inclinations between 22 and 36 deg, so this regime matters.

## 2. The constrained estimator (gaiaorb/rvincl.py v0.1.0)

The radial velocities measure the reflex amplitude s = a0 sin i (from K, P, e and the DR3 parallax) and
the star's argument of periastron omega. The estimator computes the posterior of the folded inclination
given the linear fit's Thiele-Innes estimate b and its 4x4 covariance C (other parameters profiled
exactly), using:
- isotropic orientation;
- a log-uniform true mass;
- s_obs ~ N(a0 sin i, sigma_s), with sigma_s / s = sqrt((sigma_K/K)^2 + (sigma_P/P)^2 + (e sigma_e/(1-e^2))^2);
- omega ~ N(omega_RV, sigma_omega), integrated on a Gaussian-weighted grid (spacing <= 1.5 deg, +-3 sigma);
- two omega sign conventions, omega and -omega, with equal weight. Codes differ on this; omega + 180 deg
  and a mirrored sky are already covered by the full i and Omega ranges.

The a0 integral is done in closed form. The Omega integral uses a periodic grid, refined locally for
narrow peaks and finished with Laplace's method for peaks narrower still. Inclination uses a coarse
grid, then a fine grid over the support. True mass = m sin i / sin i_f.

This is the standard joint posterior that full Bayesian codes sample: kima (Baycroft, Faria & Delisle
2026, arXiv:2606.24132), ExoPlanetLab (Hobbs 2026, arXiv:2611.13770) and orvara (Brandt et al. 2021). It is
reimplemented here, standard library only, as a deterministic computation on the frozen linear fit, so
the pre-registered pipeline does not change.

## 3. Calibration gates (epocha2_job.py, oracle mode, P <= 11 yr)

The sample is the same 300 Epochalypse systems and shard as v1. "Measured" uses the unchanged v1
verdict, so the sample does not depend on the new estimator. Oracle mode gives the true s and omega and
assumes sigma_s/s = 0.02, sigma_omega = 0. Perturbed mode applies the v1 ephemeris errors plus s and
m sin i scaled together by (1 + 0.03 N) and omega + 3 deg N, with sigma_s/s = 0.03 and sigma_omega = 3 deg.

- **C1:** constrained 95% inclination coverage >= 85% over all measured systems (needs >= 30).
- **C2:** the same coverage >= 80% for true folded inclination < 25 deg (needs >= 15).
- **C3:** constrained 95% true-mass coverage >= 85% (needs >= 30).

A constrained fit that fails to finish counts as a miss. **Outcome:** CE_PASS if all three pass,
CE_FAIL:<gates> otherwise, CE_UNDERPOWERED if any sample is too small.

Reported without gates:
- both estimators in five inclination bins;
- the perturbed mode;
- periods up to 20 yr;
- median signed and absolute errors;
- which omega branch is preferred;
- a check that the shape-only oracle coverage reproduces v1 (0.904977...) exactly.

## 4. What changes on DR4 day

Run order:
1. `rv_job.py`, then `rv_job_a1.py`, unchanged.
2. `python ce_dr4.py --points-dir /out --out-dir /out --ce-outcome <outcome> [--validated]`, on the normal
   points that rv_job_a1.py wrote.

- **If CE_PASS:** the constrained posterior becomes the primary reported inclination and true mass for
  GJ 3512 b/c, GJ 317 b/c, GJ 463 b and GJ 9773 b. It also becomes primary for the AGREES/DISAGREES checks
  of Amendment 1, section 5: a published inclination AGREES if it lies inside the 95% interval. The
  shape-only v0.2.0 numbers are reported alongside.
- **If CE_FAIL or CE_UNDERPOWERED:** the shape-only numbers stay primary. Any planet whose shape-only
  median folded inclination is below 30 deg is labelled "face-on regime: shape-only calibrated coverage
  14/21". The constrained numbers are reported as secondary.
- **GJ 9773 b:** both the printed omega (99.0) and the alternative (81.0) are run. A constrained check
  counts only if both agree, else INCONCLUSIVE_OMEGA_AMBIGUOUS (as in Amendment 1).
- **Unchanged:** the detection, consistency and true-mass verdict rules, the GJ 3512 primary run, the
  mutual inclinations (shape-only TI draws), the `--validated` gate and the scoop log.

## 5. Limits

- The 25 deg split, and so gate C2, was chosen after seeing the v1 results. C2 is evaluated on the same 21
  systems that motivated the change. It tests whether the fix works where the problem was found; it is
  not independent evidence. The perturbed mode and the bins are reported to show behaviour beyond them.
- The simulator has white Gaussian noise and one shard. Real DR4 data have CCD-level structure.
- The constrained estimator trusts the RV amplitude. A wrong K (for example from an unmodelled planet)
  biases it; the v1 consistency check (z <= 3) still guards against gross disagreement.

## 6. Credits

Yahalomi et al. 2026 (arXiv:2611.06915; Zenodo 10.5281/zenodo.23789041, CC BY 4.0) for the Epochalypse
catalogue. Baycroft, Faria & Delisle (kima) and Hobbs (ExoPlanetLab) for joint epoch-astrometry and RV posteriors,
and for showing how Thiele-Innes priors interact with inclination. Morales et al. 2025, Sozzetti 2023 and
Piccinini et al. 2026 for the planets and the predictions being tested.
