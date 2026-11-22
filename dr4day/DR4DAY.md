# Gaia DR4 release day: one job for all four near-miss tests

Fixed 2026-11-22, before Gaia DR4 (archive reopens 2026-12-02 12:00 CET).

`dr4day.py` runs, in this order and in one container:

| step | script (repository) | what it decides |
|---|---|---|
| 0 preflight | `SELECT TOP 1 source_id FROM gaiadr4.gaia_source` | stop with DR4_NOT_AVAILABLE if DR4 is not live |
| 1 validate | `validate_job.py --skip-rehearsal` (proxima-c-dr4) | VALIDATED or FAILED on the DR4 pre-release (Gaia BH3, Gaia-4, HD 114762, quasars, injections) |
| 2 proxima | `proxima_job.py` (proxima-c-dr4) | Proxima c: detected / excluded / inconclusive |
| 3 radio | `radio_job.py` (radio-planets-dr4) | TVLM 513-46546 b and GJ 896 Ab |
| 4 pceb | `etv_job.py` (pceb-planets-dr4) | QS Vir, DE CVn, RR Cae, HW Vir (2), V893 Sco |
| 5 gj3512 | `rv_job.py` (gj3512-giants-dr4) | GJ 3512 b and c true masses |
| 6 gj3512a1 | `rv_job_a1.py` (Amendment 1) | GJ 317 b/c, GJ 463 b, GJ 9773 b; DR3 predictions |
| 7 gj3512ce | `ce_dr4.py --ce-outcome X` (Amendment 2) | constrained estimator primary only if X = CE_PASS |

Rules carried over unchanged from the pre-registrations:

- `--validated` reaches steps 2-7 only if step 1's `validation.json` says `VALIDATED`.
- `X` is the `outcome` in `results_ce.json` of the `rvdr4-epochalypse-ce` calibration. If that
  calibration has not run, X is `NOT_RUN`, and the shape-only estimator stays primary.
- A failed step is logged, and the next step still runs.

Outputs: `/out/<step>/...` (each script's own files), `/out/logs/<step>.txt` and
`/out/dr4day_summary.json`.

Before filing the job, check DPAC's own DR4 orbit tables (`nss_two_body_orbit` and any exoplanet
table) for every target and log what they say (scoop log). DPAC orbits are compared with ours,
never used as inputs.

Made with AI agents working in the name of MV2 of Munim, Inc. (mv2@munimtech.com). Apache-2.0.
