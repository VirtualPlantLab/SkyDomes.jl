using SkyDomes
using Test

let

    # Check that zenith angles are calculated correctly
    n = 1000
    lats = .-π ./ 2 .+ π .* rand(n)
    DOYs = rand(1:365, n)
    decs = SkyDomes.declination.(DOYs)
    ts = 24.0 .* rand(n)
    temp = [SkyDomes.solar_angles(lat = lats[i], dec = decs[i], t = ts[i])
            for i in 1:n]
    cos_theta, theta, phi = Tuple(getindex.(temp, i) for i in 1:3)
    @test maximum(theta) <= π
    @test minimum(theta) >= 0
    @test maximum(cos_theta) <= 1.0
    @test minimum(cos_theta) >= -1.0
    @test maximum(phi) <= 2π
    @test minimum(phi) >= 0

    # Compute air mass for the different cos_theta and theta from above
    thetas = 0.0:0.01:(π / 2)
    ams = SkyDomes.air_mass.(cos.(thetas), thetas)

    # Check that extraterrestrial radiation is computed correctly
    DOYs = 1:365
    I0s = [SkyDomes.extraterrestrial(d) for d in DOYs]
    @test minimum(I0s) > 1322
    @test maximum(I0s) < 1413

    # Calculate solar radiation and components for a clear sky
    DOY = 182
    lat = π / 4
    df = 0.01
    temp = [clear_sky(lat = lat, DOY = DOY, f = x) for x in df/2:df:1-df/2]
    Igs, Idirs, Idifs = Tuple(getindex.(temp, i) for i in 1:3)
    @test minimum(Igs) > eps(Float64)
    @test all(abs.(Igs .- Idirs) .>= eps(Float64))
    @test all(abs.(Igs .- Idifs) .>= eps(Float64))
    @test all(Igs .≈ Idirs .+ Idifs)

    # Calculate daily solar radiation asusming 20% and 50% reduction with respect to clear sky
    dec = declination(DOY) # declination angle of the sun
    DL = day_length(lat, dec) # daylength in h
    Igd_sun = sum(Igs)*DL*df*3600 # Daily Ig in J/m2
    Igd_20 = Igd_sun*0.8
    Igd_50 = Igd_sun*0.5
    Iday_20 = daily_radiation(lat = lat, DOY = DOY, Igd = Igd_20)
    @test Iday_20.Igd == Igd_20
    @test Iday_20.Idif + Iday_20.Idir == Iday_20.Igd
    Iday_50 = daily_radiation(lat = lat, DOY = DOY, Igd = Igd_50)
    @test Iday_50.Igd == Igd_50
    @test Iday_20.Iod == Iday_50.Iod
    @test Iday_50.Idif + Iday_50.Idir == Iday_50.Igd

    # Calculate solar radiation and components for a cloudy sky (assume 20% reduction in Igd)
    temp = [cloudy_sky(Iday = Iday_20, lat = lat, DOY = DOY, f = x) for x in df/2:df:1-df/2]
    Igsca, Idirsc, Idifsc = Tuple(getindex.(temp, i) for i in 1:3)
    @test abs(sum(Igsca)*DL*df*3600 - Igd_20)/Igd_20 < 1e-4
    @test minimum(Igsca) > eps(Float64)
    @test all(abs.(Igsca .- Idirsc) .>= eps(Float64))
    @test all(abs.(Igsca .- Idifsc) .>= zero(Float64))
    @test all(Igsca .≈ Idirsc .+ Idifsc)
    @test all((Idifsc .- Idifs) .>= -0.5)

    # Calculate solar radiation and components for a cloudy sky (assume 50% reduction in Igd)
    temp = [cloudy_sky(Iday = Iday_50, lat = lat, DOY = DOY, f = x) for x in df/2:df:1-df/2]
    Igsca, Idirsc, Idifsc = Tuple(getindex.(temp, i) for i in 1:3)
    @test abs(sum(Igsca)*DL*df*3600 - Igd_50)/Igd_50 < 1e-4
    @test minimum(Igsca) > eps(Float64)
    @test all(abs.(Igsca .- Idirsc) .>= eps(Float64))
    @test all(abs.(Igsca .- Idifsc) .>= zero(Float64))
    @test all(Igsca .≈ Idirsc .+ Idifsc)
    @test all((Idifsc .- Idifs) .>= -4.7)

    # Calculate solar radiation and components for a cloudy sky (assume 20% reduction in Ig)
    Igsc = Igs.*0.8
    x = df/2:df:1-df/2
    temp = [cloudy_sky(Ig = Igsc[i], lat = lat, DOY = DOY, f = x[i]) for i in eachindex(x)]
    Igsca, Idirsc, Idifsc = Tuple(getindex.(temp, i) for i in 1:3)
    @test all(Igsca .== Igsc)
    @test minimum(Igsc) > eps(Float64)
    @test all(abs.(Igsc .- Idirsc) .>= eps(Float64))
    @test all(abs.(Igsc .- Idifsc) .>= zero(Float64))
    @test all(Igsc .≈ Idirsc .+ Idifsc)
    @test all((Idifsc .- Idifs) .>= -0.5) # clear sky model can lead to slightly higher dif
    @test all((Idirsc .- Idirs) .<= 0.0)

    # Calculate solar radiation and components for a cloudy sky (assume 50% reduction in Ig)
    Igsc = Igs.*0.5
    x = df/2:df:1-df/2
    temp = [cloudy_sky(Ig = Igsc[i], lat = lat, DOY = DOY, f = x[i]) for i in eachindex(x)]
    Igsca, Idirsc, Idifsc = Tuple(getindex.(temp, i) for i in 1:3)
    @test all(Igsca .== Igsc)
    @test minimum(Igsc) > eps(Float64)
    @test all(abs.(Igsc .- Idirsc) .>= eps(Float64))
    @test all(abs.(Igsc .- Idifsc) .>= zero(Float64))
    @test all(Igsc .≈ Idirsc .+ Idifsc)
    @test all((Idifsc .- Idifs) .>= -4.8) # clear sky model can lead to slightly higher dif
    @test all((Idirsc .- Idirs) .<= 0.0)

    # Test waveband conversion coefficients
    for Itype in (:direct, :diffuse)
        for waveband in (:PAR, :NIR, :UV, :red, :green, :blue)
            for mode in (:power, :flux)
                f = SkyDomes.waveband_conversion(Itype = Itype,
                    waveband = waveband,
                    mode = mode)
                @test mode == :power ? 0.0 < f < 1.0 : 0.0 < f < 6.0
            end
        end
    end

    # Check that zenith angles are calculated correctly
    lat = 0.0*pi/180
    DOY = 182
    dec = SkyDomes.declination(DOY)
    ts = 0.1:0.1:24.0
    temp = [SkyDomes.solar_angles(lat = lat, dec = dec, t = t)
            for t in ts]
    cos_thetas, thetas, phis = Tuple(getindex.(temp, i) for i in 1:3)
    day = cos_thetas .> 0
    # plot(ts[day], theta[day].*180.0./pi)
    # plot!(ts[day], phis[day].*180.0./pi)

end
