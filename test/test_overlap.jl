using SkyDomes
import PlantGraphs as PG
import PlantGeomPrimitives as PGP
import PlantGeomTurtle as PGT
import PlantRayTracer as PRT
#import GLMakie
#import PlantViz as PV
using Test
import ColorTypes: RGB
import StaticArrays: @SVector

# Simple graph that creates tiles
module Tiles
    using PlantGraphs
    using PlantRayTracer
    struct Tile <: Node
        length::Float64
        mat::Vector{Sensor{1}}
    end
end
import .Tiles

#let
    # Power and absorbed radiance are extracted with a query
    function get_power(graph, N, M)
        tiles = PG.apply(graph, PG.Query(Tiles.Tile))
        power = sum(sum(mat.power[1] for mat in tile.mat) for tile in tiles)
        area = dx*dy
        return power, power ./ area
    end

    # Auxilliary function for relative differences
    rel(x, y) = abs(x - y) / x

    ##################### Common elements ####################

    # Length of unit tile
    L = 0.5
    dx = 0.9L
    dy = L

    # Construct the tile with N triangles
    function PGT.feed!(turtle::PGT.Turtle, tile::Tiles.Tile, data)
        N = length(tile.mat)
        for i in 1:N
            PGP.Rectangle!(turtle, length = tile.length / N, width = tile.length,
                materials = tile.mat[i], colors = rand(RGB), move = true)
        end
        return nothing
    end

    # Common settings for the ray tracer (force overlap of 0.1)
    settings = PRT.RTSettings(pkill = 0.3, maxiter = 1, nx = 3, ny = 3, dx = dx, dy = dy)
    psettings = PRT.RTSettings(pkill = 0.3, maxiter = 1, nx = 3, ny = 3, parallel = true, dx = dx, dy = dy)

    # Create the ray tracing scene and run it
    function ray_trace!(scene, settings, acceleration = PRT.Naive; radiosity = 0.0,
                        nrays = 1_000_000)
        if acceleration == PRT.Naive
            acc_mesh = PRT.accelerate(scene, settings = settings, acceleration = acceleration)
        else
            acc_mesh = PRT.accelerate(scene, settings = settings, acceleration = acceleration,
                                  rule = PRT.SAH{1}(2, 5))
        end
        source = SkyDomes.sky(acc_mesh, Idir = 0.0, Idif = radiosity, nrays_dif = nrays,
                               nrays_dir = 0, sky_model = StandardSky,
                               dome_method = equal_solid_angles, ntheta = 12, nphi = 12)
        rtobj = PRT.RayTracer(acc_mesh, source)
        nrays = PRT.trace!(rtobj)
        return nothing
    end

    # Test results
    function test_results(irradiance, target = radiosity)
        @test rel(irradiance, target) < 1e-2
        return nothing
    end

    N = 10
    nw = 1
    rad = 0.5
    nrays = 100_000
    make_tile() = Tiles.Tile(L, [PRT.Sensor(1) for _ in 1:N])
    axiom = PGT.RA(-90.0) + make_tile()
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #PV.render(scene)
    minimum(getindex.(scene.vertices, 1))
    maximum(getindex.(scene.vertices, 1))
    minimum(getindex.(scene.vertices, 2))
    maximum(getindex.(scene.vertices, 2))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = rad, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(rad))

    # Check individual tiles
    function get_powers(graph)
        tiles = PG.apply(graph, PG.Query(Tiles.Tile))[1]
        powers = [mat.power[1] for mat in tiles.mat]
    end
    powers = get_powers(graph)
    @test powers[1] == 0.0
    sum(powers[2:N]) == pow
#end
