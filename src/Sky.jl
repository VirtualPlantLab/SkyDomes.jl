module Sky

import Base: getindex, lastindex, length, iterate

export UniformSky, StandardSky, CIE, radiosity,
       equal_solid_angles, equal_angle_intervals, sky, clear_sky, 
       waveband_conversion, day_length, declination

using DataFrames
using HCubature
import StaticArrays as SA
import StaticArrays: SVector
import VPL: DirectionalSource, Vec, O
import VPL.RT: AABB

include("SkyDome.jl")
include("SolarIrradiance.jl")

end
