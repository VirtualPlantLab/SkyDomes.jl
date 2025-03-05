### This file contains public API ###
# waveband_conversion
# clear_sky

#=
This code is based on:
- report SAND2012-2389 on "Global Horizontal Irradiance Clear Sky Models: Implementation and Analysis".
    - To compute solar angles
- Ineichen and Perez (2002)
    - To compute direct and diffuse solar radiation throughout the day
=#

# Calculate the declination angle
declination(DOY) = 23.45 * π / 180.0 * sin(2π / 365 * (DOY - 81))

# Calculate the solar zenith and azimuth angle and its cosine at a given location and time.
function solar_angles(; lat, t, dec)
    h = (t - 12.0) * 15.0 * π / 180.0 # hour angle of the sun
    cos_theta = cos(lat) * cos(dec) * cos(h) + sin(lat) * sin(dec) # cosine of zenith angle
    theta = acos(cos_theta) # zenith angle
    # cosine of azimuth angle for each half of the day
    cos_phi = (sin(dec)*cos(lat) - cos(h)*cos(dec)*sin(lat))/sin(theta)
    phi = acos(clamp(cos_phi, -1, 1)) # azimuth angle without correction
    phi = h < 0 ? phi : 2π - phi # compass azimuth angle with correction
    return (cos_theta, theta, phi) # return both values as a tuple
end

#=
 The extraterretrial radiation in W/m2 according to the American Society of
Civil Engineers (ASCE).
=#
function extraterrestrial(DOY)
    1367.7 * (1.0 + 0.033 * cos(DOY / 365 * 2 * π))
end

# Kasten and Young (1989) -> Note that theta is transformed to degrees
function air_mass(cos_theta, theta)
    1.0 / (cos_theta + 0.50572 * (96.07995 - theta * 180 / pi)^(-1.6354))
end

# Calculate the length of the diurnal part of a day based on latitude and declination angle
function day_length(lat, dec)
    # Calculate sunset angle with respect to solar noon, including corrections for extreme latitudes
    cosSunset = -tan(lat) * tan(dec)
    sunset = ifelse(cosSunset > 1.0, 0.0, ifelse(cosSunset < -1.0, π, acos(cosSunset)))
    # Day length in hours knowing that 1 hour = 15 degrees
    2.0 * sunset * 180.0 / π / 15.0
end

# Returns time of day in solar time (hours) as a function of day length and the
# fraction of daylength (0 = sunrise, 1 = sunset)
function timeOfDay(f, DL)
    # Time of sunrise
    tsr = 12.0 - DL / 2.0
    # diurnal time -> time of the day
    tod = tsr + f * DL
    tod
end

"""
    clear_sky(;lat, DOY, t, altitude = 0.0, TL = 4.0)

Calculate global, direct and diffuse solar radiation on the horizontal plane using
the clear sky model by Ineichen and Perez (2002).

# Arguments
- `lat`: latitude in radians
- `DOY`: day of year
- `f`: fraction of the day (0 = sunrise, 1 = sunset)
- `altitude`: altitude above sea level in meters (default 0.0)
- `TL`: Linke turbidity coefficient (default 4.0)

# Returns
A named tuple with fields:
- `Ig`: global solar radiation on the horizontal plane in W/m^2
- `Idir`: direct solar radiation on the horizontal plane in W/m^2
- `Idif`: diffuse solar radiation on the horizontal plane in W/m^2
- `theta`: solar zenith angle in radians
- `phi`: solar azimuth angle in radians

# References
Ineichen P., Perez R., A new airmass independent formulation for
the Linke turbidity coefficient, Solar Energy, Vol 73(3), pp.151–157, 2002.
"""
function clear_sky(; lat, DOY, f, altitude = 0.0, TL = 4.0)
    # Check validity of inputs
    @assert abs(lat) <= pi / 2
    @assert 0.0 <= f <= 1.0
    @assert 0 < DOY <= 365
    # Basic astronomical quantities
    Io = extraterrestrial(DOY)
    dec = declination(DOY) # declination angle of the sun
    DL = day_length(lat, dec)
    t = timeOfDay(f, DL)
    cos_theta, theta, phi = solar_angles(; lat = lat, dec = dec, t = t)
    # Clear sky model by Ineichen and Perez (2002)
    am = air_mass(cos_theta, theta)
    fh1 = exp(-altitude / 8000)
    fh2 = exp(-altitude / 1250)
    a1 = 5.09e-5 * altitude + 0.868
    a2 = 3.92e-5 * altitude + 0.0387
    # Global solar radiation on the horizontal plane
    Ig = Io * cos_theta * a1 * exp(-a2 * am * (fh1 + fh2 * (TL - 1)))
    # Direct solar radiation on the horizontal plane
    b = 0.664 + 0.163 / fh1
    Idir = Io * cos_theta * b * exp(-0.09 * am * (TL - 1))
    return (Ig = Ig, Idir = Idir, Idif = Ig - Idir, theta = theta, phi = phi)
end



"""
    cloudy_sky(Ig; lat, DOY, f)

Calculate global, direct and diffuse solar radiation on the horizontal plane using
the cloudy sky model by Spitters et al (1986) for instantaneous (not daily) measurements.

# Details

The equations to compute the fraction of incoming radiation that is diffuse were originally
developed for hourly solar radiation levels.

# Arguments
- `Ig`: Observed global solar radiation on the horizontal plane in W/m^2
- `lat`: latitude in radians
- `DOY`: day of year
- `f`: fraction of the day (0 = sunrise, 1 = sunset)

# Returns
A named tuple with fields:
- `Ig`: global solar radiation on the horizontal plane in W/m^2
- `Idir`: direct solar radiation on the horizontal plane in W/m^2
- `Idif`: diffuse solar radiation on the horizontal plane in W/m^2
- `theta`: solar zenith angle in radians
- `phi`: solar azimuth angle in radians

# References
Spitters CJ, Toussaint HA, Goudriaan J. Separating the diffuse and direct component of
global radiation and its implications for modeling canopy photosynthesis Part I. Components
of incoming radiation. Agricultural and Forest Meteorology Vol 38(1-3), pp. 217-29, 1986.
"""
function cloudy_sky(Ig; lat, DOY, f)
    # Check validity of inputs
    @assert abs(lat) <= pi / 2
    @assert 0.0 <= f <= 1.0
    @assert 0 < DOY <= 365
    @assert Ig > 0.0
    # Basic astronomical quantities
    Io = extraterrestrial(DOY)
    dec = declination(DOY) # declination angle of the sun
    DL = day_length(lat, dec)
    t = timeOfDay(f, DL)
    # Solar elevation angle
    _, theta, phi = solar_angles(; lat = lat, dec = dec, t = t)
    beta = π/2 - theta
    # Calculate Ig/Io
    IgIo =  ifelse(beta == 0.0, 0.0, min(1.0, Ig/Io/sin(beta)))
    # Compute Idf/Ig according to equation in Appendix at Spitters et al. (1986)
    R = 0.847 - 1.61*sin(beta) + 1.04*sin(beta)^2
    K = (1.47 - R)/1.66
    if IgIo <= 0.22
        IdfIg = 1.0
    elseif IgIo <= 0.35
        IdfIg = 1.0 - 6.4*(IgIo - 0.22)^2
    elseif IgIo <= K
        IdfIg = 1.47 - 1.66*IgIo
    else
        IdfIg = R
    end
    Idif = Ig*IdfIg
    Idir = Ig - Idif
    return (Ig = Ig, Idir = Idir, Idif = Idif, theta = theta, phi = phi)
end

# Convert solar irradiance to a particular waveband based on a reference solar
# spectrum for direct and diffuse solar using the Bird model. The calculations
# are included in chapter 4 of the PhD Thesis by Alejandro Morales.
# https://research.wur.nl/en/publications/dynamic-photosynthesis-under-a-fluctuating-environment-a-modellin
# Input is always in W/m2 as generated by clear_sky or similar function
const wavebands_dict = Dict(:direct => Dict(:power => Dict(:UV => 0.03,
            :blue => 0.12,
            :green => 0.16,
            :red => 0.13,
            :PAR => 0.40,
            :NIR => 0.56),
        :flux => Dict(:UV => 0.03 / 0.327,
            :blue => 0.12 / 0.263,
            :green => 0.16 / 0.217,
            :red => 0.13 / 0.184,
            :PAR => 0.40 / 0.216,
            :NIR => 0.56 / 0.105)),
    :diffuse => Dict(:power => Dict(:UV => 0.13,
            :blue => 0.25,
            :green => 0.22,
            :red => 0.13,
            :PAR => 0.60,
            :NIR => 0.27),
        :flux => Dict(:UV => 0.13 / 0.332,
            :blue => 0.25 / 0.265,
            :green => 0.22 / 0.219,
            :red => 0.13 / 0.185,
            :PAR => 0.60 / 0.226,
            :NIR => 0.27 / 0.122)))

"""
    waveband_conversion(;Itype = :direct, waveband = :PAR, mode = :power)

Returns the conversion coefficient from solar radiation (W/m2) to a give
waveband in either power (W/m2) or photon flux (umol/m2/s). The coefficients
are based on the Bird spectral model for a clear sky using June 21th in The
Netherlands (latitude 52° N).

# Arguments
- `Itype`: The type of solar radiation, either `:direct` or `:diffuse`.
- `waveband`: The waveband of interest, one of `:PAR`, `:UV`, `:blue`, `:red`, `:green`, or `:NIR`.
- `mode`: The physical units of the target, either `:power` (W/m^2) or `:flux` (umol/m^2/s).

# Examples
```julia
waveband_conversion(Itype = :diffuse, waveband = :UV, mode = :flux)
waveband_conversion(waveband = :NIR)
```
"""
function waveband_conversion(; Itype = :direct, waveband = :PAR, mode = :power)
    # Check the inputs for validity
    @assert Itype in (:direct, :diffuse)
    @assert mode in (:power, :flux)
    @assert waveband in (:PAR, :UV, :blue, :red, :green, :NIR)
    # Extract the coversion coefficient from the nested dictionary
    wavebands_dict[Itype][mode][waveband]
end
