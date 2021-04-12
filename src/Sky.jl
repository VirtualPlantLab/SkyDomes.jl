module Sky

import Base: getindex, lastindex, length, iterate

export UniformSky, StandardSky, CIE, radiance, radiosity, polar,
       equal_solid_angles, equal_angle_intervals

using DataFrames
using RCall 
using HCubature
import Unitful
const UN = Unitful
import Unitful: °, rad, hr, ustrip, @u_str, m, mbar, minute, d, s, K
import StaticArrays
const SA = StaticArrays
import StaticArrays: SVector

include("Utils.jl")
include("SkyDome.jl")
include("SolarAngles/SPA.jl")

end
