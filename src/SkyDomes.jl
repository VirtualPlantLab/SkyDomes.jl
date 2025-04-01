module SkyDomes

import Base: getindex, lastindex, length, iterate

export UniformSky, StandardSky, CIE, radiosity,
    equal_solid_angles, equal_angle_intervals, sky, clear_sky, cloudy_sky,
    waveband_conversion, day_length, declination, daily_radiation

using DataFrames
using HCubature
import StaticArrays as SA
import StaticArrays: SVector
import PlantRayTracer: DirectionalSource, AccMesh, AABB, center
import PlantGeomPrimitives: Mesh, Vec

include("Sources.jl")
include("SolarIrradiance.jl")

end
