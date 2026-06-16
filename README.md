# SkyDomes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://virtualplantlab.com/stable/SkyDomes/API/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://virtualplantlab.com/dev/SkyDomes/API/)
[![CI](https://github.com/VirtualPlantLab/SkyDomes.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/VirtualPlantLab/SkyDomes.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/VirtualPlantLab/SkyDomes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/VirtualPlantLab/SkyDomes.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/605131958.svg)](https://zenodo.org/doi/10.5281/zenodo.10256565)

The package SkyDomes provides a function to calculate the solar radiation on a
horizontal surface (for clear skies) as a function of latitude, day of year and
time of the day and for different wavebands. In addition, it can generate light
sources as required by the [Virtual Plant Lab](https://github.com/VirtualPlantLab/VirtualPlantLab.jl) to
simulate the light distribution in a 3D scene.

## Installation

To install SkyDomes.jl, you can use the following command:

```julia
] add SkyDomes
```

Or, if you prefer the development version:

```julia
import Pkg
Pkg.add(url = "https://github.com/VirtualPlantLab/SkyDomes.jl.git", rev = "master")
```

## Usage


### Solar radiation

Use the `clear_sky` function to calculate the solar radiation on a horizontal
plane as a function of day of year, latitude (in degrees) and the relative solar
time of the day (`f = 0` is sunrise, `f = 1` is sunset). The function returns
a named tuple with the total, direct and diffuse solar radiation in W/m² as well
as the solar zenith and azimuth angles in degrees. For example:

```julia
using SkyDomes
lat = 52.0 # latitude in degrees
DOY = 182
f = 0.5 # solar noon
Ig, Idir, Idif, theta, phi = clear_sky(lat = lat, DOY = DOY, f = f) 
```

The values `Ig`, `Idir` and `Idif` are the total, direct and diffuse solar
radiation in W/m². The angles `theta` and `phi` are the solar zenith and azimuth
angles in degrees; they are needed later to point the direct light source in the
right direction. The function `waveband_conversion` can be used to convert
broadband irradiance to specific wavebands (UV, PAR, NIR, blue, green or red) as
well as converting from W/m² to µmol/m²/s, assuming particular spectra for
direct and diffuse solar radiation. For example:

```julia
f_PAR_dir = waveband_conversion(Itype = :direct, waveband = :PAR, mode = :flux)
Idir_PAR = f_PAR_dir*Idir # PAR direct in µmol/m²/s
f_PAR_dif = waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :flux)
Idif_PAR = f_PAR_dif*Idif # PAR diffuse in µmol/m²/s
```

### Light sources for ray tracing

Once the direct and diffuse solar radiation in the relevant wavebands and units
have been calculated, the function `sky` can be used to generate the light
sources required by VPL to simulate the light distribution in a 3D scene. For
example, a simple horizontal tile (representing soil) in VPL may be created as
follows:

```julia
using VirtualPlantLab
import GLMakie
r = Rectangle(length = 2.0, width = 1.0)
rotatey!(r, 90.0) # To put it in the XY plane
VirtualPlantLab.translate!(r, Vec(0.0, 0.5, 0.0))
```

To use the ray tracer 3D scene requires adding optical properties (a `Lambertian` material here) and
linking them to the mesh using `add!`. We keep a reference to the material object
so we can read the absorbed power after tracing. FOr 3D visualization we just need a color object.
For complicated scenes see the
[VPL documentation](http://virtualplantlab.com/) for examples:

```julia
import ColorTypes: RGB
soil = Lambertian(τ = 0.0, ρ = 0.2)
scene = Mesh()
add!(scene, r, materials = soil, colors = RGB(0.8, 0.7, 0.5))
render(scene)
```

The `sky` function requires a scene that has been processed by `accelerate` in
order to set up the grid cloner (see the
[grid cloner guide](../../howto/GridCloner.md) for details). The clone spacing
should match the scene dimensions:

```julia
settings = RTSettings(parallel = true, nx = 3, ny = 3, dx = 2.0, dy = 1.0)
acc_scene = accelerate(scene, settings = settings)
```

If we want to compute the amount of solar radiation absorbed by the tile, we
need to create a series of light sources. The function `sky` is used for
that purpose. Note that the solar zenith (`theta`) and azimuth (`phi`) angles
returned by `clear_sky` are passed directly to orient the direct light source:

```julia
sources = sky(acc_scene,
             Idir = Idir_PAR,          # Direct solar radiation from above
             nrays_dir = 1_000_000,    # Number of rays for direct solar radiation
             theta_dir = theta,        # Solar zenith angle (degrees)
             phi_dir = phi,            # Solar azimuth angle (degrees)
             Idif = Idif_PAR,          # Diffuse solar radiation from above
             nrays_dif = 10_000_000,   # Total number of rays for diffuse solar radiation
             sky_model = StandardSky,  # Angular distribution of solar radiation
             dome_method = equal_solid_angles, # Discretization of the sky dome
             ntheta = 9,               # Number of discretization steps in the zenith angle
             nphi = 12)                # Number of discretization steps in the azimuth angle
```

The function takes the accelerated scene as input to ensure that light sources
scale with the scene. Direct solar radiation is represented by a single
directional light source that will emit a number of rays given by `nrays_dir`.
Diffuse solar radiation is represented by a hemispherical dome of directional
light sources that will emit a total of `nrays_dif` rays. The angular
distribution of the diffuse solar radiation and the discretization of the sky dome
can be modified via `dome_method`, `sky_model`, `ntheta` and `nphi`. See API
documentation and [VPL documentation](http://virtualplantlab.com/) for details.

Once the light sources are created, a ray tracing object can be generated
combining all the elements above:

```julia
rtobj = RayTracer(acc_scene, sources, settings = settings);
```

And the ray tracing can be performed by calling the `trace!` function:

```julia
trace!(rtobj);
```

As expected, the amount of solar radiation absorbed by the tile equals the
total in the scene:

```julia
power(soil)[1]/area(r) ≈ (Idir_PAR + Idif_PAR)*0.8
```

See [VPL documentation](http://virtualplantlab.com/) for more details and
tutorials on ray tracing simulations.
