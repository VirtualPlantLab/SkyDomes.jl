using Sky
using Test

let

    # Check that zenith angles are calculated correctly    
    n = 1000
    lats = .-π ./ 2 .+ π .* rand(n)
    DOYs = rand(1:365, n)
    decs = Sky.declination.(DOYs)
    ts = 24.0 .* rand(n)
    temp = [Sky.solar_zenith_angle(lat = lats[i], dec = decs[i], t = ts[i]) for i in 1:n]
    cos_theta, theta = Tuple(getindex.(temp, i) for i in 1:2)
    @test maximum(theta) <= π
    @test minimum(theta) >= 0
    @test maximum(cos_theta) <= 1.0
    @test minimum(cos_theta) >= -1.0

    # Compute air mass for the different cos_theta and theta from above
    thetas = 0.0:0.01:(π / 2)
    ams = Sky.air_mass.(cos.(thetas), thetas)

    # Check that extraterrestrial radiation is computed correctly
    DOYs = 1:365
    I0s = [Sky.extraterrestrial(d) for d in DOYs]
    @test minimum(I0s) > 1322
    @test maximum(I0s) < 1413

    # Calculate solar radiation and components for a clear sky
    temp = [Sky.clear_sky.(lat = π / 4, DOY = 182, f = x) for x in 0.01:0.01:0.99]
    Igs, Idirs, Idifs = Tuple(getindex.(temp, i) for i in 1:3)
    @test minimum(Igs) > eps(Float64)
    @test all(abs.(Igs .- Idirs) .>= eps(Float64))
    @test all(abs.(Igs .- Idifs) .>= eps(Float64))
    @test all(Igs .≈ Idirs .+ Idifs)

    # Test waveband conversion coefficients
    for Itype in (:direct, :diffuse)
        for waveband in (:PAR, :NIR, :UV, :red, :green, :blue)
            for mode in (:power, :flux)
                f = Sky.waveband_conversion(Itype = Itype, waveband = waveband, mode = mode)
                @test mode == :power ? 0.0 < f < 1.0 : 0.0 < f < 6.0
            end
        end
    end
end
