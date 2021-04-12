using Unitful
import Sky
using Plots
using StaticArrays
import Unitful: rad, °, mbar, m, °C, K, d, s

# Example in manual
year  = 2003
month = 10
day   = 17
hour  = 12
minute = 30
second = 30
timezone = -7
longitude = -105.1786°
latitude =  39.742476°
elevation = 1830.14m 
pressure = 820mbar
temperature = convert(typeof(1.0K), 11°C)
slope = 30°
azm_rotation = -10°


ΔT = 67s
Sky.sunrise_sunset_transit(year, month, day, latitude, longitude, ΔT).*24
test = Sky.SPA(year, month, day, hour, minute, second, timezone, latitude, longitude, 
          slope, azm_rotation,elevation, pressure, temperature, ΔT)

# Difference between my implementation and reference          
ϵ_ang = convert(typeof(1.0rad), 0.0001°)
ϵ_time = 0.5/3600/24
abs(test.θ - 50.11162°) < ϵ_ang # Zenith
abs(test.Φ - 194.34024°) < ϵ_ang # Azimuth
abs(test.I - 25.18700°) < ϵ_ang # Incidence
abs(test.R - (13/24 + 12/60/24 + 43.46/3600/24)) < ϵ_time # Sunrise
abs(test.T - (18/24 + 46/60/24 + 4.97/3600/24)) < ϵ_time # Transit
abs(test.S - (20/60/24 + 19.19/3600/24)) < ϵ_time # Sunset

@code_warntype Sky.SPA(year, month, day, hour, minute, second, timezone, latitude, longitude, 
slope, azm_rotation,elevation, pressure, temperature, ΔT)

@code_warntype Sky.geocentric_coordinates(year, month, 1.0, ΔT)

using BenchmarkTools

Sky.heliocentric_latitude(1.0)
@benchmark Sky.heliocentric_radius(1.0)



ccall((:deg2rad, "C:/juliabin/spa.o"), Float64, (Float64, ), 90.0)
