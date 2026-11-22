# gj3512-giants-dr4

**How heavy are the two giant planets around one of the smallest stars known to host them?**

GJ 3512 is a red dwarf with an eighth of the Sun's mass, 9.5 parsecs away. Radial velocities
show two gas giants around it (Morales et al. 2019, *Science*; Morales et al. 2025, A&A,
arXiv:2507.15516). Core accretion has trouble forming giants around so small a star, and only
their minimum masses (0.46 and 0.45 Jupiter masses) are known. The discoverers wrote that Gaia's
epoch astrometry in Data Release 4 would settle the true masses.

This repository fixes, **before DR4 is released on 2026-12-02**, exactly how that measurement
will be made and graded:

- Rules and verdicts: [PREREGISTRATION.md](PREREGISTRATION.md)
- Predictions: the star must wobble by at least **0.127 mas** (planet b, 203 d) and
  **0.636 mas** (planet c, 6.4 yr); the true wobble is these divided by sin i.
- Code: `gaiaorb/rvdr4.py` (an exact joint linear fit of both RV orbits to Gaia's along-scan
  epochs), `rv_job.py` (the DR4 run), `test_rvdr4.py` (tests on synthetic epochs).
- The shared Gaia DR4 code (`gaiaorb/`, `io_gaia.py`, `etv_job.py`) is the same as in
  [proxima-c-dr4](https://github.com/munimtechnologies/proxima-c-dr4),
  [radio-planets-dr4](https://github.com/munimtechnologies/radio-planets-dr4) and
  [pceb-planets-dr4](https://github.com/munimtechnologies/pceb-planets-dr4).
- Run after release: `python -m unittest -v test_rvdr4 test_gaiaorb` then
  `python rv_job.py --data-dir /data --out-dir /out [--validated]`

Status (2026-11-22): pre-registered; waits for Gaia DR4.

## Credits

Morales et al. 2019 (Science 365, 1441) and Morales et al. 2025 (arXiv:2507.15516) for the
planets, the orbital elements and the forecast; Lopez-Santiago et al. 2020 (AJ 160, 273) for the
first evidence of planet c; ESA's Gaia mission and DPAC for the astrometry.

## Disclosure and licence

Made with AI agents working in the name of MV2 of Munim, Inc. (mv2@munimtech.com). No human
expert has checked it yet. Code: Apache-2.0 (see LICENSE). Text and results: CC BY 4.0.
