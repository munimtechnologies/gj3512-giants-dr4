# Scoop and prior-work check, Amendment 1 (2026-11-22)

Searched the web, arXiv and the Gaia archive on 2026-11-22 for GJ 317, GJ 463 and GJ 9773 true
masses, astrometric orbits and Gaia DR4.

- **Gaia DR4** is not out (release 2026-12-02, 12:00 CET). No DR4 orbit can exist yet. Gaia DR3
  has no non-single-star solution for any of the three stars (non_single_star = 0).
- **GJ 317:**
  - Anglada-Escude et al. 2012 (ApJ 746, 37; arXiv:1111.2623): ground-based astrometry; b's
    mass < 2x minimum at 99%, inclination > 25 deg.
  - Morales et al. 2025: DR3 excess noise, assuming coplanarity; i = 33.9 (+6.1/-1.4) deg,
    3.03 and 2.56 M_Jup.
  - Nothing else found from Hipparcos-Gaia or epoch astrometry.
- **GJ 463 b:**
  - Sozzetti 2023 (A&A 670, L17; arXiv:2302.00413): PMa + RV; i = 27 +- 3 / 152 (+2/-3) deg,
    3.4 +- 0.3 / 3.6 +- 0.4 M_Jup. Says DR4 should give the full orbit.
  - Piccinini et al. 2026 (A&A, arXiv:2601.09401): used as a benchmark; i = 27.7 / 152.3 deg,
    3.3 / 3.4 M_Jup.
  - Morales et al. 2025: 21.8 deg, 4.44 M_Jup.
  - These disagree by about 2 sigma. No DR4 result exists.
- **GJ 9773 b:** discovered by Morales et al. 2025 (PMa mass 1.43 +0.52/-0.08 M_Jup at 22.3 deg);
  no other analysis found.
- **Other DR4-preparation projects:** coracle-us/gaia-dr4 on GitHub (outer companions of hot and
  warm Jupiter hosts) does not overlap. Yahalomi et al. 2026 (Epochalypse) is a simulation
  benchmark that reports no orbit-fitter results; we use it for calibration.
- **Noted in passing:** Morales et al. 2025 Table 4 prints e cos(omega) = +0.034, e sin(omega) =
  0.214 and omega = 99.0 deg for GJ 9773 b, which are mutually inconsistent (they imply omega =
  81.0 deg). Both are run; see the amendment, section 2.

Re-check on release day, before publishing any result: the DR4 NSS / exoplanet tables, arXiv
astro-ph.EP listings of 2-4 December 2026, and the authors' pages.
