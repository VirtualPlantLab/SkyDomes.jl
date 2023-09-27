using Sky
using Test
using Documenter
import Aqua

# Test examples on documentation (jldoctest blocks)
DocMeta.setdocmeta!(Sky,
    :DocTestSetup,
    :(using Sky);
    recursive = true)
doctest(Sky)

# Aqua
@testset "Aqua" begin
    Aqua.test_all(Sky, ambiguities = false, project_extras = false)
    Aqua.test_ambiguities([Sky])
end

@testset "Solar irradiance" begin
    include("test_solar_irradiance.jl")
end
@testset "Irradiance distribution" begin
    include("test_irradiance_distribution.jl")
end
@testset "Ray tracing" begin
    include("test_raytracing.jl")
end
