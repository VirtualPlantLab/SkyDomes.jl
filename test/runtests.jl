using SkyDomes
using Test
using Documenter
import Aqua

# Test examples on documentation (jldoctest blocks)
DocMeta.setdocmeta!(SkyDomes,
    :DocTestSetup,
    :(using SkyDomes);
    recursive = true)
doctest(SkyDomes)

# Aqua
@testset "Aqua" begin
    Aqua.test_all(SkyDomes, ambiguities = false, project_extras = false)
    Aqua.test_ambiguities([SkyDomes])
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
@testset "Ray tracing tiles" begin
    include("test_raytracer_tiles.jl")
end
