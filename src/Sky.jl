module Sky

import Base: getindex, lastindex, length, iterate

export UniformSky, StandardSky, CIE, radiance, radiosity,
       equal_solid_angles, equal_angle_intervals, sky, clear_sky, 
       waveband_conversion

using DataFrames
using HCubature
import Unitful
const UN = Unitful
import Unitful: °, rad, hr, ustrip, @u_str, m, mbar, minute, d, s, K
import StaticArrays
const SA = StaticArrays
import StaticArrays: SVector
import VPL: DirectionalSource, Vec, O
import VPL.RT: AABB

include("Utils.jl")
include("SkyDome.jl")
include("SolarIrradiance.jl")

end
