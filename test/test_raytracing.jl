using Sky
using PlantGeomPrimitives
using PlantRayTracer
using Test

let

    # Steps to integrate the sky into VPL
    # 0. Create the scene in VPL
    r = Rectangle(length = 2.0, width = 1.0)
    rotatey!(r, -π/2) # To put it in the XY plane
    translate!(r, Vec(0.0, 0.5, 0.0))
    r2 = deepcopy(r)
    translate!(r2, Vec(0.0, 0.0, -1.0))
    materials = [Black(2), Black(2)]
    ids = [1, 1, 2, 2]
    scene = Scene(mesh = Mesh([r,r2]), material_ids = ids, materials = materials)

    # 1. Generate the solar irradiance at a specific moment in time
    lat = 52.0*π/180.0
    DOY = 182
    f = 0.5
    Ig, Idir, Idif = Sky.clear_sky(lat = lat, DOY = DOY, f = f)

    # 2. Convert to different wavebands (PAR, NIR, UV, red, green, blue)
    Idir_PAR = Sky.waveband_conversion(Itype = :direct, waveband = :PAR, mode = :flux)*Idir
    Idif_PAR = Sky.waveband_conversion(Itype = :diffuse, waveband = :PAR, mode = :flux)*Idif

    # 3. Create a sky dome of directional light sources using the solar irradiance
    sources = Sky.sky(scene, Idir = (Idir_PAR, Idir_PAR), Idif = (Idif_PAR, Idif_PAR),
                      nrays_dif = 10_000_000,  nrays_dir = 1_000_000,
                      sky_model = StandardSky,
                      dome_method = equal_solid_angles, ntheta = 9, nphi = 12);

    # 4. Create a ray tracing object with the sky dome as the light sources
    settings = RTSettings(nx = 3, ny = 3, dx = 2.0, dy = 1.0, parallel = true)
    rtobj = RayTracer(scene, sources, settings = settings);

    # 5. Run the ray tracer!
    nrays_traced = trace!(rtobj)

    # 6. Check the power in the scene
    power_out = sum(source.power*source.nrays for source in sources)[1]
    @test materials[1].power[1]/power_out ≈ 1.0
    @test materials[1].power[2]/power_out ≈ 1.0
    @test materials[2].power[1]/power_out ≈ 0.0
    @test materials[2].power[2]/power_out ≈ 0.0
    @test materials[1].power[1]/area(r) ≈ Idir_PAR + Idif_PAR
    @test materials[1].power[2]/area(r) ≈ Idir_PAR + Idif_PAR

end
