# SkyDomes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://virtualplantlab.com/stable/SkyDomes/API/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://virtualplantlab.com/dev/SkyDomes/API/)
[![CI](https://github.com/VirtualPlantLab/SkyDomes.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/VirtualPlantLab/SkyDomes.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/VirtualPlantLab/SkyDomes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/VirtualPlantLab/SkyDomes.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

The package SkyDomes provides a function to calculate the solar radiation on a
horizontal surface (for clear skies) as a function of latitude, day of year and
time of the day and for different wavebands. In addition, it can generate light
sources as required by [VirtualPlantLab.jl](https://github.com/VirtualPlantLab/VirtualPlantLab.jl) to
simulate the light distribution in a 3D scene.

## Installation

To install SkyDomes.jl, you can use the following command:

```julia
] add SkyDomes
```

Or, if you prefer the development version:

```julia
import Pkg
Pkg.add(url = "https://github.com/VirtualPlantLab/SkyDomes.jl", rev = "master")
```

## Usage


### Solar radiation

Use the `clear_sky` function to calculate the solar radiation on a horizontal
plane as a function of day of year, latitude (in radians) and the relative solar
time of the day (`f = 0` is sunrise, `f = 1` is sunset). The function returns
the total solar radiation in W/m² as well as direct and diffuse components. For
example:

```julia
using SkyDomes
lat = 52.0*π/180.0 # latitude in radians
DOY = 182
f = 0.5 # solar noon
Ig, Idir, Idif, theta, phi = clear_sky(lat = lat, DOY = DOY, f = f) # W/m2
```

The values `Ig`, `Idir` and `Idif` are the total, direct and diffuse solar
radiation in W/m². The function `waveband_conversion` can be used to convert
these values to specific wavebands (UV, PAR, NIR, blue, green or red) as well
as converting from W/m² to umol/m²/s, assuming particular spectra for
direct and diffuse solar radiation. For example:

```julia
f_PAR_dir = waveband_conversion(Itype = :direct, waveband = :PAR, mode = :flux)
Idir_PAR = f_PAR_dir*Idir # PAR in umol/m²/s
f_PAR_dif = waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :flux)
Idif_PAR = f_PAR_dif*Idif # PAR in umol/m²/s
```

### Light sources for ray tracing

Once the direct and diffuse solar radiation in the relevant wavebands and units
have been calculated, the function `SkyDomes` can be used to generate the light
sources required by VPL to simulate the light distribution in a 3D scene. For
example, a simple horizontal tile (representing soil) in VPL may be created as
follows:

```julia
using VirtualPlantLab
r = Rectangle(length = 0.5, width = 0.5)
rotatey!(r, -π/2) # To put it in the XY plane
translate!(r, Vec(0.0, 0.5, 0.0))
import GLMakie
render(r)
```

A 3D scene requires adding optical properties (e.g., a black material property)
and linking these to the mesh (for complicated scenes see [VPL documentation](http://virtualplantlab.com/)
for examples):

```julia
materials = [Black()]
scene = Scene(mesh = r, material_ids = [1,1], materials = materials)
```

If we want to compute the amount of solar radiation absorbed by this tile, we
need to create a series of light sources. The function `SkyDomes` can be used for
that purpose:

```julia
sources = sky(scene,
             Idir = Idir_PAR, # Direct solar radiation from above
             nrays_dir = 1_000_000, # Number of rays for direct solar radiation
             theta_dir = theta, # Solar zenith angle for direct solar radiation
             phi_dir = phi, # Solar azimuth angle for direct solar radiation
             Idif = Idif_PAR, # Diffuse solar radiation from above
             nrays_dif = 10_000_000, # Total number of rays for diffuse solar radiation
             sky_model = StandardSky, # Angular distribution of solar radiation
             dome_method = equal_solid_angles, # Discretization of the SkyDomes dome
             ntheta = 9, # Number of discretization steps in the zenith angle
             nphi = 12) # Number of discretization steps in the azimuth angle
render!(sources) # To visualize the dome of light sources
```

The function takes the scene as input to ensure that light sources scale with
the scene. Direct solar radiation is represented by a single directional light
source that will emit a number of rays given by `nrays_dir`. Diffuse solar
radiation is represented by a hemispherical dome of directional light sources
that will emit a total of `nrays_dif` rays. The angular distribution of the
diffuse solar radiation and the discretization of the SkyDomes dome can be modified
via `dome_method`, `sky_model`, `ntheta` and `nphi`. See API documentation and
[VPL documentation](http://virtualplantlab.com/) for details.

Once the light sources are created, a ray tracing object can be generated
combining all the elements above:

```julia
settings = RTSettings(parallel = true)
rtobj = RayTracer(scene, sources, settings = settings);
```

And the ray tracing can be performed by calling the `trace!` function:

```julia
trace!(rtobj);
```

As expected, the amount of solar radiation absorbed by the tile equals the
total in the scene (`Ig`):

```julia
power(materials[1])[1]/area(r) ≈ Idir_PAR + Idif_PAR
```

See [VPL documentation](http://virtualplantlab.com/) for more details and
tutorials on ray tracing simulations

## Roadmap for future development

The package is still under development. The following features are planned:

- Calculate fraction of direct and diffuse radiation from measurements of actual
solar radiation.

- Calculate the solar radiation on a tilted surface.

- Allow for a different orientation of the 3D scene from VPL.

- More advanced algorithms for solar radiation and spectrum.

See Issues for additional features to be implemented.
