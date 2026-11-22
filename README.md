# gj3512-giants-dr4

**How heavy are the giant planets around small red dwarfs, and were the Gaia DR3 estimates right?**

GJ 3512 is a red dwarf with an eighth of the Sun's mass, 9.5 parsecs away. Radial velocities
show two gas giants around it (Morales et al. 2019, *Science*; Morales et al. 2025, A&A,
arXiv:2507.15516). Core accretion has trouble forming giants around so small a star, and only
their minimum masses (0.46 and 0.45 Jupiter masses) are known. The discoverers wrote that Gaia's
epoch astrometry in Data Release 4 would settle the true masses.

This repository fixes, **before DR4 is released on 2026-12-02**, exactly how those measurements
will be made and graded:

- **GJ 3512 b and c** (v0.1.0): [PREREGISTRATION.md](PREREGISTRATION.md). The star must wobble by
  at least **0.127 mas** (b, 203 d) and **0.636 mas** (c, 6.4 yr); the true wobble is these
  divided by sin i.
- **Amendment 1, GJ 317 b and c, GJ 463 b, GJ 9773 b** (v0.2.0):
  [PREREGISTRATION-amendment1.md](PREREGISTRATION-amendment1.md). These are the other M-dwarf
  giants of Morales et al. 2025. The amendment also tests their published Gaia DR3-based
  inclinations, checks whether the two-planet systems are coplanar as assumed, and decides between
  two disagreeing GJ 463 b solutions: about 22 deg and 4.4 M_Jup (Morales et al. 2025) versus
  about 27-28 deg and 3.3-3.6 M_Jup (Sozzetti 2023; Piccinini et al. 2026).
- **Calibration v1, graded 2026-11-22: EPOCHALYPSE_PASS.**
  [VALIDATION-epochalypse.md](VALIDATION-epochalypse.md) set the rules and
  [GRADING-epochalypse-v1.md](GRADING-epochalypse-v1.md) reports the result. The frozen fitter was run
  on 300 systems of the Epochalypse simulated DR4 epoch astrometry (Yahalomi et al. 2026,
  arXiv:2611.06915).
  - False alarms: 0 of 300.
  - Inclination intervals: 90.5% cover the truth.
  - Inconsistent with the RVs: 2.2%.
  - A diagnostic done afterwards found the weak spot: nearly face-on orbits. Below 25 deg only 14 of 21
    intervals covered the truth, with inclinations biased toward edge-on by a median 11 deg.
- **Amendment 2** (2026-11-22, before DR4): [PREREGISTRATION-amendment2.md](PREREGISTRATION-amendment2.md).
  It adds an inclination estimator that also uses the RV amplitude a0 sin i (`gaiaorb/rvincl.py`),
  calibrates it on the same systems under pre-stated gates (`epocha2_job.py`), and fixes the rule that
  decides which estimator is primary on DR4 day (`ce_dr4.py`).
- **I/O fix** (2026-11-22, before DR4): [notes/io-fix-2026-11-22.md](notes/io-fix-2026-11-22.md).
  The Gaia Archive renamed its DataLink data structures on 2026-10-14 (INDIVIDUAL became
  DATAMODEL_STANDARD). `io_gaia.py` 0.2.0 asks for the new name first and rejects error responses.
  No analysis changes.
- **DR4 day:** [dr4day/DR4DAY.md](dr4day/DR4DAY.md). `dr4day.py` runs all four of MV2's DR4 near-miss
  tests in one job, in a fixed order: validation, Proxima c, the radio planets, the eclipse-timing
  planets, then this repository's three steps.

Code:
- `gaiaorb/rvdr4.py`: an exact joint linear fit of the RV orbits to Gaia's along-scan epochs.
- `gaiaorb/rvincl.py`: the RV-amplitude-constrained inclination posterior (Amendment 2).
- Runs:
  - `rv_job.py` (GJ 3512) and `rv_job_a1.py` (Amendment 1);
  - `ce_dr4.py` (Amendment 2);
  - `epocha_job.py` and `epocha2_job.py` (calibrations);
  - `dr4day.py` (all four DR4 tests in one job).
- Tests on synthetic epochs: `test_rvdr4.py`, `test_rvdr4_a1.py`, `test_epocha.py`, `test_rvincl.py`,
  `test_epocha2.py`, `test_ce_dr4.py`, `test_io_gaia.py`, `test_dr4day.py`.

The shared Gaia DR4 code (`gaiaorb/`, `io_gaia.py`, `etv_job.py`) is the same as in
[proxima-c-dr4](https://github.com/munimtechnologies/proxima-c-dr4),
[radio-planets-dr4](https://github.com/munimtechnologies/radio-planets-dr4) and
[pceb-planets-dr4](https://github.com/munimtechnologies/pceb-planets-dr4).

To run after release (one job): `python dr4day.py --data-dir /data --out-dir /out --ce-outcome <outcome>`.
The same steps can also be run one at a time:
1. `python -m unittest -v test_rvdr4 test_rvdr4_a1 test_rvincl test_ce_dr4 test_gaiaorb test_io_gaia`
2. `python rv_job.py --data-dir /data --out-dir /out [--validated]`
3. `python rv_job_a1.py --data-dir /data --out-dir /out [--validated]`
4. `python ce_dr4.py --points-dir /out --out-dir /out --ce-outcome <outcome of epocha2_job.py> [--validated]`

Status (2026-11-22): pre-registered. Calibration v1 passed. The Amendment 2 calibration and a live
DataLink probe are queued. The measurements wait for Gaia DR4.

## Credits

- Morales et al. 2019 (Science 365, 1441) and Morales et al. 2025 (arXiv:2507.15516): the planets,
  orbital elements, DR3-based masses and the DR4 forecast.
- Lopez-Santiago et al. 2020 (AJ 160, 273): the first evidence of GJ 3512 c.
- Johnson et al. 2007 and Anglada-Escude et al. 2012 (ApJ 746, 37): GJ 317 b and its astrometry.
- Endl et al. 2022: GJ 463 b.
- Sozzetti 2023 (A&A 670, L17) and Piccinini et al. 2026 (arXiv:2601.09401): the GJ 463 b
  dynamical masses.
- Yahalomi et al. 2026 (arXiv:2611.06915; Zenodo 10.5281/zenodo.23789041, CC BY 4.0): the
  Epochalypse catalogue.
- Baycroft, Faria & Delisle 2026 (kima, arXiv:2606.24132), Hobbs 2026 (ExoPlanetLab,
  arXiv:2611.13770) and Brandt et al. 2021 (orvara): joint epoch-astrometry and RV posteriors, on
  which Amendment 2's estimator is modelled.
- ESA's Gaia mission and DPAC: the astrometry.

## Disclosure and licence

Made with AI agents working in the name of MV2 of Munim, Inc. (mv2@munimtech.com). No human
expert has checked it yet. Code: Apache-2.0 (see LICENSE). Text and results: CC BY 4.0.
