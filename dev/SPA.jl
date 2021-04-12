
# Wrapper around the SPA program
module SPA

# Struct passed to spa_calculate
mutable struct spa_data
    # Inputs
    year::Cint
    month::Cint
    day::Cint
    hour::Cint
    minute::Cint
    second::Cdouble
    delta_ut1::Cdouble
    delta_t::Cdouble
    timezone::Cdouble
    longitude::Cdouble
    latitude::Cdouble
    elevation::Cdouble
    pressure::Cdouble
    temperature::Cdouble
    slope::Cdouble
    azm_rotation::Cdouble
    atmos_refract::Cdouble
    fun::Cint
    # Intermediate values that are filled in during the calculations
    jd::Cdouble
    jc::Cdouble
    jde::Cdouble
    jce::Cdouble
    jme::Cdouble
    l::Cdouble
    b::Cdouble
    r::Cdouble
    theta::Cdouble
    beta::Cdouble
    x0::Cdouble
    x1::Cdouble
    x2::Cdouble
    x3::Cdouble
    x4::Cdouble
    del_psi::Cdouble
    del_epsilon::Cdouble
    epsilon0::Cdouble
    epsilon::Cdouble
    del_tau::Cdouble
    lamda::Cdouble
    nu0::Cdouble
    nu::Cdouble
    alpha::Cdouble
    delta::Cdouble
    h::Cdouble
    xi::Cdouble
    del_alpha::Cdouble
    delta_prime::Cdouble
    alpha_prime::Cdouble
    h_prime::Cdouble
    e0::Cdouble
    del_e::Cdouble
    e::Cdouble
    eot::Cdouble
    srha::Cdouble
    ssha::Cdouble
    sta::Cdouble
    # Final output
    zenith::Cdouble
    azimuth_astro::Cdouble
    azimuth::Cdouble
    incidence::Cdouble
    suntransit::Cdouble
    sunrise::Cdouble
    sunset::Cdouble
end

# Facilitate creation of struct
create_spa_data(;year = 2003, month = 10, day = 17, hour = 12, minute = 30, 
                 second = 30, delta_ut1 = 0.0, delta_t = 67.0, 
                timezone = -7.0, longitude = -105.1786, latitude = 39.742476, 
                elevation = 1830.14, pressure = 820.0, temperature = 11.0, 
                slope = 30.0, azm_rotation = -10.0, atmos_refract = 0.5667, 
                fun = 3) =

    spa_data(year, month, day, hour, minute, second, delta_ut1, delta_t, 
             timezone, longitude, latitude, elevation, pressure, 
             temperature, slope, azm_rotation, atmos_refract, fun,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# Wrappers to all the functions in SPA

# double deg2rad(double degrees);
deg2rad(degrees::Cdouble) = ccall((:deg2rad, "C:/juliabin/spa.o"), Cdouble, (Cdouble, ), degrees)

# double rad2deg(double radians);
rad2deg(radians::Cdouble) = ccall((:rad2deg, "C:/juliabin/spa.o"), Cdouble, (Cdouble, ), radians)

# double limit_degrees(double degrees);
limit_degrees(degrees::Cdouble) = ccall((:limit_degrees, "C:/juliabin/spa.o"), Cdouble, (Cdouble, ), degrees)

# double third_order_polynomial(double a, double b, double c, double d, double x);
third_order_polynomial(a::Cdouble, b::Cdouble, c::Cdouble, d::Cdouble, x::Cdouble) = 
    ccall((:third_order_polynomial, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), 
           a, b, c, d, x)

# double geocentric_right_ascension(double lamda, double epsilon, double beta);
geocentric_right_ascension(lamda::Cdouble, epsilon::Cdouble, beta::Cdouble) = 
    ccall((:geocentric_right_ascension, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble), lamda, epsilon, beta)

# double geocentric_declination(double beta, double epsilon, double lamda);
geocentric_declination(beta::Cdouble, epsilon::Cdouble, lamda::Cdouble) = 
    ccall((:geocentric_declination, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble), beta, epsilon, lamda)

# double observer_hour_angle(double nu, double longitude, double alpha_deg);
observer_hour_angle(nu::Cdouble, longitude::Cdouble, alpha_deg::Cdouble) = 
    ccall((:observer_hour_angle, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble), nu, longitude, alpha_deg)

# void   right_ascension_parallax_and_topocentric_dec(double latitude, double elevation,
# 	         double xi, double h, double delta, double *delta_alpha, double *delta_prime);
right_ascension_parallax_and_topocentric_dec(latitude::Cdouble, elevation::Cdouble, 
                            xi::Cdouble, h::Cdouble, delta::Cdouble, 
                            delta_alpha::Vector{Cdouble}, delta_prime::Vector{Cdouble}) = 
    ccall((:right_ascension_parallax_and_topocentric_dec, "C:/juliabin/spa.o"), 
    Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}), 
    latitude, elevation, xi, d, h, delta, delta_alpha, delta_prime)

# double topocentric_right_ascension(double alpha_deg, double delta_alpha);
topocentric_right_ascension(alpha_deg::Cdouble, delta_alpha::Cdouble) = 
    ccall((:topocentric_right_ascension, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble), alpha_deg, delta_alpha)

# double topocentric_local_hour_angle(double h, double delta_alpha);
topocentric_local_hour_angle(h::Cdouble, delta_alpha::Cdouble) = 
    ccall((:topocentric_local_hour_angle, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble), h, delta_alpha)

# double topocentric_elevation_angle(double latitude, double delta_prime, double h_prime);
topocentric_elevation_angle(latitude::Cdouble, delta_alpha::Cdouble, h_prime::Cdouble) = 
    ccall((:topocentric_elevation_angle, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble), latitude, delta_alpha, h_prime)

# double atmospheric_refraction_correction(double pressure, double temperature,
#                                          double atmos_refract, double e0);
atmospheric_refraction_correction(pressure::Cdouble, temperature::Cdouble, 
                                  atmos_refract::Cdouble, e0::Cdouble) = 
    ccall((:atmospheric_refraction_correction, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble, Cdouble, Cdouble), 
          pressure, temperature, atmos_refract, e0)

# double topocentric_elevation_angle_corrected(double e0, double delta_e);
topocentric_elevation_angle_corrected(e0::Cdouble, delta_e::Cdouble) = 
    ccall((:topocentric_elevation_angle_corrected, "C:/juliabin/spa.o"), Cdouble, 
          (Cdouble, Cdouble), e0, delta_e)

# double topocentric_zenith_angle(double e);
topocentric_zenith_angle(e::Cdouble) = 
    ccall((:topocentric_zenith_angle, "C:/juliabin/spa.o"), Cdouble, (Cdouble,), e)

# double topocentric_azimuth_angle_astro(double h_prime, double latitude, double delta_prime);
topocentric_azimuth_angle_astro(h_prime::Cdouble, latitude::Cdouble, delta_prime::Cdouble) = 
    ccall((:topocentric_azimuth_angle_astro, "C:/juliabin/spa.o"), Cdouble, 
    (Cdouble, Cdouble, Cdouble), h_prime, latitude, delta_prime)

# double topocentric_azimuth_angle(double azimuth_astro);
topocentric_azimuth_angle(azimuth_astro::Cdouble) = 
    ccall((:topocentric_azimuth_angle, "C:/juliabin/spa.o"), Cdouble, 
    (Cdouble, ), azimuth_astro)

# int spa_calculate(spa_data *spa);
spa_calculate(data::spa_data) = ccall((:spa_calculate, "C:/juliabin/spa.o"), 
                                      Cint, (Ref{spa_data}, ), [data])

end

import .SPA

# Generate 1000 different situations by sampling within the hypercube of 
# parameters and put it in a tuple with the corresponding names
function add_names(values)
    (year = Int(floor(values[1])), month = Int(floor(values[2])), 
    day = Int(floor(values[3])), hour = Int(floor(values[4])), 
    minute = Int(floor(values[5])), second = values[6], delta_ut1 = values[7], 
    delta_t = values[8], timezone = values[9], longitude = values[10], 
    latitude = values[11], elevation = values[12], pressure = values[13], 
    temperature = values[14], slope = values[15], azm_rotation = values[16])
end

using Sobol
lb = [1500, 1, 0, 0, 0, 0, -1, -8000, -18, -180, -90, -1000, 600, -50, -360, -360]
ub = [2500, 12, 31, 23, 60, 60, 1, 8000, 18, 180, 90, 7000, 1300, 50, 360, 360]
qrng = SobolSeq(lb, ub)
skip(qrng, 1024)

scenarios = Vector{SPA.spa_data}(undef, 1000)
for i in 1:length(scenarios)
    while(true)
        sample = next!(qrng)
        sample = add_names(sample)
        scenario = SPA.create_spa_data(;sample...)
        flag = SPA.spa_calculate(scenario)
        if flag == 0 
            scenarios[i] = scenario
            break
        end
    end
end

# Save the file for usage in package tests
using BSON

bson("./test/SPA.bson", Dict(:scenarios => scenarios))

# Compare with my own implementation
import Sky
using Unitful
import Unitful: °, m, mbar, L, °C, K, s
using Test

SPA_scenarios = BSON.load("./test/SPA.bson")[:scenarios]

rel(a,b) = abs((a - b)/b)

c = 0
for scenario in SPA_scenarios
   global c += 1
#    scenario = SPA_scenarios[12]
   println(c)
    julia_ver = Sky.SPA(scenario.year, scenario.month, scenario.day, 
                          scenario.hour, scenario.minute, scenario.second, 
                          scenario.timezone, scenario.latitude*°, scenario.longitude*°, 
                          scenario.slope*°, scenario.azm_rotation*°, 
                          scenario.elevation*m, scenario.pressure*mbar, 
                          convert(typeof(1.0K), scenario.temperature*°C), 
                          scenario.delta_t*s)


    @test abs(convert(typeof(1.0°), julia_ver.θ) - scenario.zenith*°) < 0.005° #2eps(Float64)^0.25
    @test abs(convert(typeof(1.0°), julia_ver.Φ) - scenario.azimuth*°) < 0.005° #2eps(Float64)^0.25
    @test abs(convert(typeof(1.0°), julia_ver.Γ) - scenario.azimuth_astro*°) < 0.005° #2eps(Float64)^0.25
    @test abs(convert(typeof(1.0°), julia_ver.I) - scenario.incidence*°) < 0.005° #2eps(Float64)^0.25
    @test rel(julia_ver.R, scenario.sunrise) < 2eps(Float64)^0.25
    @test rel(julia_ver.S, scenario.sunset) < 2eps(Float64)^0.25
    @test rel(julia_ver.T, scenario.suntransit) < 2eps(Float64)^0.25
end

import Unitful: rad, °

1rad > 1°
-0.010525204964702635 > -(0.26667 + Rₛᵣ)°
# Test the individual functions
import .SPA
import Sky
using Sobol
using Test

# geocentric_right_ascension
sob = SobolSeq([-180,23.2,-90.0], [180.0,23.6,90.0])
for i in 1:100
    angles = next!(sob)
    lamda = angles[1]
    epsilon = angles[2]
    beta = angles[3]
    rlamda = lamda*π/180
    repsilon = epsilon*π/180
    rbeta = beta*π/180
    C = SPA.geocentric_right_ascension(lamda, epsilon, beta) 
    @show C
    Julia = Sky.geogentric_sun_ascension(rlamda, repsilon, rbeta)*180/π
    @test C ≈ Julia
end


# geocentric_declination
sob = SobolSeq([-90,23.2,-180.0], [90.0,23.6,180.0])
for i in 1:100
    angles = next!(sob)
    lamda = angles[1]
    epsilon = angles[2]
    beta = angles[3]
    rlamda = lamda*π/180
    repsilon = epsilon*π/180
    rbeta = beta*π/180
    C = SPA.geocentric_declination(beta, epsilon, lamda) 
    Julia = Sky.geocentric_sun_declination(rlamda, repsilon, rbeta)*180/π
    @test C ≈ Julia
end

# observer_hour_angle
sob = SobolSeq([0,-180.0, 0.0], [360,180.0, 360.0])
for i in 1:100
    angles = next!(sob)
    nu = angles[1]
    longitude = angles[2]
    alpha = angles[3]
    rnu = nu*π/180
    rlongitude = longitude*π/180
    ralpha = alpha*π/180
    C = SPA.observer_hour_angle(nu, longitude, alpha) 
    Julia = Sky.local_hour_angle(rlongitude, rnu, ralpha)*180/π
    @test C ≈ Julia
end


