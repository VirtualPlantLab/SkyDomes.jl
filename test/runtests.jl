using Sky
using Test

@testset "Solar position" begin include("test_solar_position.jl") end
@testset "Irradiance distribution" begin include("test_irradiance_distribution.jl") end

