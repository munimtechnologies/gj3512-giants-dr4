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
- **Calibration before DR4:** [VALIDATION-epochalypse.md](VALIDATION-epochalypse.md). The frozen
  fitter is run on the Epochalypse simulated DR4 epoch astrometry (Yahalomi et al. 2026,
  arXiv:2611.06915), which contains injected companions of known orbits and zero-companion nulls.
  Three gates are fixed in advance: false-alarm rate, inclination-interval coverage and
  consistency.

Code:
- `gaiaorb/rvdr4.py`: an exact joint linear fit of the RV orbits to Gaia's along-scan epochs.
- Runs: `rv_job.py` (GJ 3512), `rv_job_a1.py` (Amendment 1), `epocha_job.py` (calibration).
- Tests on synthetic epochs: `test_rvdr4.py`, `test_rvdr4_a1.py`, `test_epocha.py`.

The shared Gaia DR4 code (`gaiaorb/`, `io_gaia.py`, `etv_job.py`) is the same as in
[proxima-c-dr4](https://github.com/munimtechnologies/proxima-c-dr4),
[radio-planets-dr4](https://github.com/munimtechnologies/radio-planets-dr4) and
[pceb-planets-dr4](https://github.com/munimtechnologies/pceb-planets-dr4).

To run after release: first `python -m unittest -v test_rvdr4 test_rvdr4_a1 test_gaiaorb`, then
`python rv_job.py --data-dir /data --out-dir /out [--validated]` and
`python rv_job_a1.py --data-dir /data --out-dir /out [--validated]`.

Status (2026-11-22): pre-registered; the calibration is queued; the measurements wait for Gaia DR4.

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
- ESA's Gaia mission and DPAC: the astrometry.

## Disclosure and licence

Made with AI agents working in the name of MV2 of Munim, Inc. (mv2@munimtech.com). No human
expert has checked it yet. Code: Apache-2.0 (see LICENSE). Text and results: CC BY 4.0.
