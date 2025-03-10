using SkyDomes
using PlantGeomPrimitives
import PlantGeomPrimitives as PGP
using PlantRayTracer
using Test

let

    # Steps to integrate the sky into VPL
    # 0. Create the mesh in VPL
    r = Rectangle(length = 2.0, width = 1.0)
    rotatey!(r, -π / 2) # To put it in the XY plane
    translate!(r, Vec(0.0, 0.5, 0.0))
    r2 = deepcopy(r)
    translate!(r2, Vec(0.0, 0.0, -1.0))
    mats = [Black(2), Black(2)]
    ids = [1, 1, 2, 2]
    mesh = Mesh([r, r2])
    PGP.add_property!(mesh, :materials, mats[ids])

    # 1. Generate the solar irradiance at a specific moment in time
    lat = 52.0 * π / 180.0
    DOY = 182
    f = 0.5
    Ig, Idir, Idif, theta = SkyDomes.clear_sky(lat = lat, DOY = DOY, f = f)

    # 2. Convert to different wavebands (PAR, NIR, UV, red, green, blue)
    Idir_PAR = SkyDomes.waveband_conversion(Itype = :direct,
        waveband = :PAR,
        mode = :flux) *
               Idir
    Idif_PAR = SkyDomes.waveband_conversion(Itype = :diffuse,
        waveband = :PAR,
        mode = :flux) *
               Idif

    # 3. Create a sky dome of directional light sources using the solar irradiance
    settings = RTSettings(nx = 10, ny = 10, dx = 2.0, dy = 1.0, parallel = true)
    acc_mesh = accelerate(mesh, settings = settings)
    sources = SkyDomes.sky(acc_mesh, Idir = (Idir_PAR, Idir_PAR), Idif = (Idif_PAR, Idif_PAR),
        nrays_dif = 10_000_000, nrays_dir = 1_000_000, theta_dir = theta,
        sky_model = StandardSky,
        dome_method = equal_solid_angles, ntheta = 9, nphi = 12)

    # 4. Create a ray tracing object with the sky dome as the light sources
    rtobj = RayTracer(acc_mesh, PlantRayTracer.materials(mesh), sources, settings)

    # 5. Run the ray tracer!
    nrays_traced = trace!(rtobj)

    # 6. Check the power in the mesh
    power_out = sum(source.power * source.nrays for source in sources)[1]
    @test power(mats[1])[1] / power_out ≈ 1.0
    @test power(mats[1])[2] / power_out ≈ 1.0
    @test power(mats[2])[1] / power_out ≈ 0.0
    @test power(mats[2])[2] / power_out ≈ 0.0
    @test power(mats[1])[1] / area(r) ≈ Idir_PAR + Idif_PAR
    @test power(mats[1])[2] / area(r) ≈ Idir_PAR + Idif_PAR

    # Modify so the cloner is smaller than the soil tile
    settings = RTSettings(nx = 10, ny = 10, dx = 1.5, dy = 0.75, parallel = true)
    acc_mesh = accelerate(mesh, settings = settings)
    sources = SkyDomes.sky(acc_mesh, Idir = (Idir_PAR, Idir_PAR), Idif = (Idif_PAR, Idif_PAR),
                            nrays_dif = 10_000_000, nrays_dir = 1_000_000, theta_dir = theta,
                            sky_model = StandardSky,
                            dome_method = equal_solid_angles, ntheta = 9, nphi = 12)
    rtobj = RayTracer(acc_mesh, PlantRayTracer.materials(mesh), sources, settings)
    nrays_traced = trace!(rtobj)
    power_out = sum(source.power * source.nrays for source in sources)[1]
    @test power(mats[1])[1] / power_out ≈ 1.0
    @test power(mats[1])[2] / power_out ≈ 1.0
    @test power(mats[2])[1] / power_out ≈ 0.0
    @test power(mats[2])[2] / power_out ≈ 0.0
    # Use the area defined by the grid cloner (1.5 x 0.75) rather than the total area
    @test power(mats[1])[1] / 1.5/0.75 ≈ Idir_PAR + Idif_PAR
    @test power(mats[1])[2] / 1.5/0.75 ≈ Idir_PAR + Idif_PAR
end
