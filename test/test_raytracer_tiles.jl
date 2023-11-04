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
struct Tile{N, M} <: Node
    length::Float64
    mat::Vector{Sensor{M}}
end
end
import .Tiles

let
    # Power and absorbed radiance are extracted with a query
    function get_power(graph, N, M)
        tiles = PG.apply(graph, PG.Query(Tiles.Tile{N, M}))
        power = sum(sum(sum(mat.power) for mat in tile.mat)
                    for tile in tiles)
        area = length(tiles) * tiles[1].length^2
        return power, power ./ area
    end

    # Auxilliary function for relative differences
    rel(x, y) = abs(x - y) / x

    ##################### Common elements ####################

    # Length of unit tile
    L = 0.5

    # Construct the tile with N triangles
    function PGT.feed!(turtle::PGT.Turtle, tile::Tiles.Tile{N, M}, data) where {N, M}
        for _ in 1:N
            PGP.Rectangle!(turtle, length = tile.length / N, width = tile.length,
                material = tile.mat[1], color = rand(RGB), move = true)
        end
        return nothing
    end

    # Common settings for the ray tracer
    settings = PRT.RTSettings(pkill = 0.3, maxiter = 1, nx = 3, ny = 3)
    psettings = PRT.RTSettings(pkill = 0.3, maxiter = 1, nx = 3, ny = 3, parallel = true)

    # Create the ray tracing scene and run it
    function ray_trace!(scene,
        settings,
        acceleration = PRT.Naive;
        radiosity = 0.0,
        nrays = 1000)
        source = SkyDomes.sky(scene, Idir = 0.0, Idif = radiosity, nrays_dif = nrays,
                               nrays_dir = 0, sky_model = StandardSky,
                               dome_method = equal_solid_angles, ntheta = 9, nphi = 12)
        if acceleration == PRT.Naive
            rtobj = PRT.RayTracer(scene,
                source,
                settings = settings,
                acceleration = acceleration)
        else
            rtobj = PRT.RayTracer(scene,
                source,
                settings = settings,
                acceleration = PRT.BVH,
                rule = PRT.SAH{1}(2, 5))
        end
        nrays = PRT.trace!(rtobj)
        return nothing
    end

    # Test results
    function test_results(irradiance, target = radiosity)
        @test rel(irradiance, target) < 1e-2
        return nothing
    end

    ##################### One tile with two triangles and single wavelength ####################

    # Construct the scene
    N = 1
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    axiom = PGT.RA(90.0) + Tiles.Tile{N, nw}(L, [PRT.Sensor(1)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Scene(graph)
    #render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    ##################### One tile with multiple triangles and single wavelength ####################
    N = 10
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    axiom = PGT.RA(90.0) + Tiles.Tile{N, nw}(L, [PRT.Sensor(1)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Scene(graph)
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle

    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    ##################### One tile with multiple triangles and multiple wavelength ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000
    axiom = PGT.RA(90.0) +
            Tiles.Tile{N, nw}(L, [PRT.Sensor(2)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Scene(graph)
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    ##################### Many tiles with multiple triangles and single wavelength ####################
    N = 10
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Sensor(1)])
    axiom = PGT.RA(90.0) + make_tile() + make_tile()
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Scene(graph)
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle

    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))


    ##################### Many tiles with multiple triangles and multiple wavelengths ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Sensor(2)])
    axiom = PGT.RA(90.0) + make_tile() + make_tile()
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Scene(graph)
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity))


    ##################### Two graphs - many tiles with multiple triangles and multiple wavelengths ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Sensor(2)])
    axiom1 = PGT.RA(90.0) + make_tile() + make_tile()
    graph1 = PG.Graph(axiom = axiom)
    axiom2 = PGT.RA(90.0) + PGT.F(2.0) + make_tile() + make_tile()
    graph2 = PG.Graph(axiom = axiom2)
    scene = PGP.Scene([graph1, graph2])
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity))


    ##################### Two graphs - many tiles with multiple triangles and multiple wavelengths - use add! ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 500_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Sensor(2)])
    axiom1 = PGT.RA(90.0) + make_tile() + make_tile()
    graph1 = PG.Graph(axiom = axiom1)
    axiom2 = PGT.RA(90.0) + PGT.F(2.0) + make_tile() + make_tile()
    graph2 = PG.Graph(axiom = axiom2)
    scene = PGP.Scene(graph1)
    scene_extra = PGP.Scene(graph2)
    mat = PRT.Sensor(2)
    PGP.add!(scene, mesh = scene_extra.mesh, material = mat, color = rand(RGB))
    #PV.render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra.mesh)
    test_results(irradiance, sum(radiosity))

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra.mesh)
    test_results(irradiance, sum(radiosity))

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra.mesh)
    test_results(irradiance, sum(radiosity))

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra.mesh)
    test_results(irradiance, sum(radiosity))

    # Source is at an angle
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra.mesh)
    test_results(irradiance, sum(radiosity))
end
