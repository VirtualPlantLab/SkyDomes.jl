using Sky
using Test

@testset "Solar irradiance" begin include("test_solar_irradiance.jl") end
@testset "Irradiance distribution" begin include("test_irradiance_distribution.jl") end
@testset "Ray tracing" begin include("test_raytracing.jl") end
