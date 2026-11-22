# Scoop and prior-art check for Amendment 2 (2026-11-22)

Searched on 2026-11-22: web search for Gaia DR4 epoch astrometry, face-on inclination bias and joint
RV fits; web search for GJ 463 b true mass and inclination; arXiv abstract pages of the papers below.

- **Joint RV + epoch-astrometry posteriors are standard.** Our estimator is not a new method.
  - kima (Baycroft, Faria & Delisle 2026, arXiv:2606.24132, A&A): samples in cos i because uniform or
    Gaussian Thiele-Innes priors imply a face-on-favouring inclination prior in joint fits. In its
    face-on test the RVs tighten the inclination bound sharply.
  - ExoPlanetLab (Hobbs 2026, arXiv:2611.13770): 180 injections; reports 90% coverage of 95% intervals
    for periods. Inclination coverage and face-on bias are not discussed. No GJ 3512, GJ 317, GJ 463 or
    GJ 9773 analysis.
  - pyGaiafits on Gaia-4b (arXiv:2611.15537) and orvara (Brandt et al. 2021) combine RVs with
    astrometry as well.
- **GJ 463 b:** still Sozzetti 2023 (152 +2/-3 deg, i.e. a folded 28 deg, 3.6 +- 0.4 M_Jup), Piccinini
  et al. 2026 (compatible) and Morales et al. 2025 (21.8 deg). No DR4-based value has been posted (DR4 is
  not public).
- **What is ours:** the calibration of one frozen, pre-registered RV-fixed linear DR4 fitter on an
  independent simulator, the measured face-on under-coverage of its shape-only inclination (14 of 21),
  and the pre-registered switch to an RV-amplitude-constrained estimator under stated gates, all before
  DR4.
