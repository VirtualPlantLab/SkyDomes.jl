# SkyDomes release notes

We started keeping track of changes in the `NEWS.md` file after version 0.1.10.

# SkyDomes 1.0.1 (2026-06-11)

All internal calculations now use degrees throughout, completing the migration
started in 1.0.0 (which only changed the public API). The note in 1.0.0 that
"internal calculations remain in radians" no longer applies.

* All trigonometric calls replaced with degree variants (`cosd`, `sind`, `tand`,
  `acosd`).
* Radian constants (`π/2`, `2π`) replaced with their degree equivalents (`90.0`,
  `360.0`) in discretisation routines and radiosity formulas.
* Normalisation denominators in `radiosity` updated from `/ π` to `/ 180` to
  match the degree-valued sector widths.
* All internal degree↔radian conversions removed from `solar_angles`,
  `air_mass`, `day_length`, `clear_sky`, `daily_radiation`, and `cloudy_sky`.
* **Exception**: inside the CIE radiance formula (`fCIE`) the angular distance
  `γₛ` is kept as a radian-valued local variable because the CIE model
  parameters (`b`, `d`) were calibrated for radian angles and appear in
  non-trigonometric terms (`exp(d·γₛ)`). All externally stored and passed
  angles remain in degrees.

# SkyDomes 1.0.0 (2026-06-10)

**Breaking changes**: All public API functions now express angles in degrees instead of
radians (both inputs and outputs). The internal calculations remain in radians.

The affected functions and their angle arguments/outputs are:

* `sky`: inputs `theta_dir`, `phi_dir`, `α`, `alpha_soil`, and `beta_soil` now in degrees.
  Default values updated accordingly (`α = 180.0`, `beta_soil = 180.0`).
* `clear_sky`: input `lat` and outputs `theta` and `phi` now in degrees.
* `daily_radiation`: input `lat` now in degrees.
* `cloudy_sky`: input `lat` and outputs `theta` and `phi` now in degrees.
* `CIE`: inputs `θₛ` and `Φₛ` now in degrees.
* `declination`: return value now in degrees (was radians).
* `day_length`: inputs `lat` and `dec` now in degrees.

# SkyDomes 0.2.0 (2026-06-10)

* Update creation of the sky dome to account to support crop rows not aligned North-South,
  and tilted, oriented slopes. All the calculations are done in PlantRayTracer.jl

# SkyDomes 0.1.10 (2026-01-14)

* Update dependencies and make sure it works on Julia 1.12
