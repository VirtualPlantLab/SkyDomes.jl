#= 
    This code is based on report NREL/TP-560-34302 on "Solar Position Algorithm 
    for Solar Radiation Applications".

    To avoid errors, users should always pass quantities to the functions rather
than plain values. Any units can be used (as long as they are compatible with
the type of variable in question). All internal calculations are type as indicated
in the manual, but conversions are not implemented in the code, but rather
handled implicitly by Unitful.

    Also, whereas hard-coded constants are implemented as Float64, time and
coordinates are not typed, in order to be compatible with, e.g., dual numbers.
=#

################################################################################
#################################### Time ######################################
################################################################################

#=
    Difference between Terrestrial Time and Universal Time due 
cumulative effect of deviations between Earth's rotation and the length of day 
in atomic time. These polynomials are taken from NASA Eclipse website
(https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html)
Note that year is time expressed in units of year in the Gregorian calendar
=#
function UT1_to_TT(y)
    if y < -500
        u = (y - 1820)/100
        return (-20 + 32u^2)s
    elseif y < 500
        u = y/100
        return (10583.6 - 1014.41u + 33.78311u^2 - 5.952053u^3 -
                 0.1798452u^4 + 0.022174192u^5 + 0.0090316521u^6)s
    elseif y < 1600
        u = (y - 1000)/100
        return (1574.2 - 556.01u + 71.23472u^2 + 0.319781u^3 -
                  0.8503463u^4 - 0.005050998u^5 + 0.0083572073u^6)s
    elseif y < 1700
        t = y - 1600
        return (120 - 0.9808t - 0.01532t^2 + t^3/7129)s
    elseif y < 1800
        t = y - 1700
        return (8.83 + 0.1603t - 0.0059285t^2 + 0.00013336t^3 - t^4/1174000)s
    elseif y < 1860
        t = y - 1800
        return (13.72 - 0.332447t + 0.0068612t^2 + 0.0041116t^3 - 0.00037436t^4 +
                 0.0000121272t^5 - 0.0000001699t^6 + 0.000000000875t^7)s
    elseif y < 1900
        t = y - 1860
        return (7.62 + 0.5737t - 0.251754t^2 + 0.01680668t^3 - 0.0004473624t^4 + 
                 t^5/233174)s
    elseif y < 1920
        t = y - 1900
        return (-2.79 + 1.494119t - 0.0598939t^2 + 0.0061966t^3 - 0.000197t^4)s
    elseif y < 1941
        t = y - 1920
        return (21.20 + 0.84493t - 0.076100t^2 + 0.0020936t^3)s
    elseif y < 1961
        t = y - 1950
        return (29.07 + 0.407t - t^2/233 + t^3/2547)s
    elseif y < 1986
        t = y - 1975
        return (45.45 + 1.067*t - t^2/260 - t^3/718)s
    elseif y < 2005
        t = y - 2000
        return (63.86 + 0.3345t - 0.060374t^2 + 0.0017275t^3 + 0.000651814t^4 +
                 0.00002373599t^5)s
    elseif y < 2050
        t = y - 2000
        return (62.92 + 0.32217t + 0.005589t^2)s
    elseif y < 2150
        return (-20 + 32*((y-1820)/100)^2 - 0.5628*(2150 - y))s
    else
        u = (y - 1820)/100
        return (-20 + 32u^2)s
    end

end

# Calculate Julian day on Gregorian calendar (Eq. 4)
function julian_day(year, month, dom, Gregorian::Bool = true)
    if month < 3
        year = year - 1
        month = month + 12
    end
    # Calculate Julian day (Eq. 4)
    JD = floor(365.25*(year + 4716)) + floor(30.6001*(month + 1)) + 
           dom - 1524.5
    if JD > 2299160
        A = year ÷ 100
        B = 2 - A + (A ÷ 4)
        JD = JD + B
    end
    return JD
end

# Julian century (Eq. 6)
JD_to_JC(JD) = (JD - 2451545)/36525

# Julian Ephemeris Day (Eq. 5)
JD_to_JDE(JD, ΔT) = JD + ΔT/86400s

# Julian Ephemeris century (Eq. 7)
JDE_to_JCE(JDE) = (JDE - 2451545)/36525

# Julian Ephemeris millenium (Eq. 8)
JCE_to_JME(JCE) = JCE/10

leap_months = cumsum(SVector(0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
no_leap_months = cumsum(SVector(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

# Calculate month and dom from fractional doy
function month_day(doy, leap::Bool = false)
    if leap
        month = findfirst(leap_months .>= doy) - 1
        dom = doy - leap_months[month]
    else
        month = findfirst(no_leap_months .>= doy) - 1
        dom = doy - no_leap_months[month]        
    end
    return month, dom
end

# Check for leap years
function is_leap(year)
    iszero(rem(year, 4)) && iszero(rem(year, 100)) && iszero(rem(year, 400)) && 
    (return true)
    return false
end

################################################################################
################################# Astronomy ####################################
################################################################################

#=
    Heliocentric coordinates: Position of Earth in spherical coordinate system 
                              centered at the Sun.

    Geocentric coordinates: Position of Sun in spherical coordinate system
                              centered at the Earth.
=#

# Auxilliary function to constrain angle to the circle (used in different locations)
function constrain_angle_circle(α)
    limit = 2π*rad
    F = α/limit
    L = limit*(F - floor(F))
    return L < zero(typeof(L)) ? L + limit : L
end

function constrain_angle_half_circle(α)
    limit = π*rad
    F = α/limit
    L = limit*(F - floor(F))
    return L < zero(typeof(L)) ? L + limit : L
end

function constrain_angle_half_circle_pm(α)
    limit = 2π*rad
    F = α/limit
    L = limit*(F - floor(F))
    if L < -limit/2
        L = L + limit
    elseif L > limit/2
        L = L - limit
    end
    return L
end

let
    # Table A4.2 Earth Periodic Terms in report
    A = (SVector(175347046, 3341656, 34894, 3497, 3418, 3136, 2676, 2343, 1324,
                 1273, 1199, 990, 902, 857, 780, 753, 505, 492, 357, 317, 284,
                 271, 243, 206, 205, 202, 156, 132, 126, 115, 103, 102, 102, 
                 99, 98, 86, 85, 85, 80, 79, 75, 74, 74, 70, 62, 61, 57, 56, 56,
                 52, 52, 51, 49, 41, 41, 39, 37, 37, 36, 36, 33, 30, 30, 25.0),
         SVector(628331966747, 206059, 4303, 425, 119, 109, 93, 72, 68, 67, 59,
                 56, 45, 36, 29, 21, 19, 19, 17, 16, 16, 15, 12, 12, 12, 12,
                 11, 10, 10, 9, 9, 8, 6, 6.0),
         SVector(52919, 8720, 309, 27, 16, 16, 10, 9, 7, 5, 4, 4, 3, 3, 3, 3, 3, 
                 3, 2, 2.0),
         SVector(289, 35, 17, 3, 1, 1, 1.0),
         SVector(114, 8, 1.0),
         SVector(1.0),
         SVector(280.0, 102, 80, 44, 32),
         SVector(9.0, 6),
         SVector(100013989.0, 1670700, 13956, 3084, 1628, 1576, 925, 542, 472,
                 346, 329, 307, 243, 212, 186, 175, 110, 98, 86, 86, 65, 63, 57, 
                 56, 49, 47, 45, 43, 39, 38, 37, 37, 36, 35, 33, 32, 32, 28, 28, 
                 26),
         SVector(103019.0, 1721, 702, 32, 31, 25, 18, 10, 9, 9),
         SVector(4359.0, 124, 12, 9, 6, 3),
         SVector(145.0, 7),
         SVector(4))

    B = (SVector(0, 4.6692568, 4.6261, 2.7441, 2.8289, 3.6277, 4.4181, 6.1352, 
                 0.7425, 2.0371, 1.1096, 5.233, 2.045, 3.508, 1.179, 2.533, 
                 4.583, 4.205, 2.92, 5.849, 1.899, 0.315, 0.345, 4.806, 1.869,
                 2.458, 0.833, 3.411, 1.083, 0.645, 0.636, 0.976, 4.267, 6.21,
                 0.68, 5.98, 1.3, 3.67, 1.81, 3.04, 1.76, 3.5, 4.68, 0.83, 
                 3.98, 1.82, 2.78, 4.39, 3.47, 0.19, 1.33, 0.28, 0.49, 5.37, 
                 2.4, 6.17, 6.04, 2.57, 1.71, 1.78, 0.59, 0.44, 2.74, 3.16),
         SVector(0, 2.678235, 2.6351, 1.59, 5.796, 2.966, 2.59, 1.14, 1.87, 
                 4.41, 2.89, 2.17, 0.4, 0.47, 2.65, 5.34, 1.85, 4.97, 2.99, 
                 0.03, 1.43, 1.21, 2.83, 3.26, 5.27, 2.08, 0.77, 1.3, 4.24, 
                 2.7, 5.64, 5.3, 2.65, 4.67),
         SVector(0, 1.0721, 0.867, 0.05, 5.19, 3.68, 0.76, 2.06, 0.83, 4.66, 
                 1.03, 3.44, 5.14, 6.05, 1.19, 6.12, 0.31, 2.28, 4.38, 3.75),
         SVector(5.844, 0, 5.49, 5.2, 4.72, 5.3, 5.97),
         SVector(3.142, 4.13, 3.84),
         SVector(3.14),
         SVector(3.199, 5.422, 3.88, 3.7, 4),
         SVector(3.9, 1.73),
         SVector(0, 3.0984635, 3.05525, 5.1985, 1.1739, 2.8469, 5.453, 4.564, 
                 3.661, 0.964, 5.9, 0.299, 4.273, 5.847, 5.022, 3.012, 5.055, 
                 0.89, 5.69, 1.27, 0.27, 0.92, 2.01, 5.24, 3.25, 2.58, 5.54, 
                 6.01, 5.36, 2.39, 0.83, 4.9, 1.67, 1.84, 0.24, 0.18, 1.78, 
                 1.21, 1.9, 4.59),
         SVector(1.10749, 1.0644, 3.142, 1.02, 2.84, 1.32, 1.42, 5.91, 1.42, 0.27),
         SVector(5.7846, 5.579, 3.14, 3.63, 1.87, 5.47),
         SVector(4.273, 3.92),
         SVector(2.56))

    C = (SVector(0, 6283.07585, 12566.1517, 5753.3849, 3.5231, 77713.7715,
                 7860.4194, 3930.2097, 11506.7698, 529.691, 1577.3435, 5884.927,
                 26.298, 398.149, 5223.694, 5507.553, 18849.228, 775.523, 
                 0.067, 11790.629, 796.298, 10977.079, 5486.778, 2544.314,
                 5573.143, 6069.777, 213.299, 2942.463, 20.775, 0.98, 
                 4694.003, 15720.839, 7.114, 2146.17, 155.42, 161000.69, 
                 6275.96, 71430.7, 17260.15, 12036.46, 5088.63, 3154.69, 
                 801.82, 9437.76, 8827.39, 7084.9, 6286.6, 14143.5, 
                 6279.55, 12139.55, 1748.02, 5856.48, 1194.45, 8429.24, 
                 19651.05, 10447.39, 10213.29, 1059.38, 2352.87, 6812.77,
                 17789.85, 83996.85, 1349.87, 4690.48),
         SVector(0, 6283.07585, 12566.1517, 3.523, 26.298, 1577.344, 18849.23,
                 529.69, 398.15, 5507.55, 5223.69, 155.42, 796.3, 775.52, 7.11,
                 0.98, 5486.78, 213.3, 6275.96, 2544.31, 2146.17, 10977.08, 
                 1748.02, 5088.63, 1194.45, 4694, 553.57, 6286.6, 1349.87,
                 242.73, 951.72, 2352.87, 9437.76, 4690.48),
         SVector(0, 6283.0758, 12566.152, 3.52, 26.3, 155.42, 18849.23, 
                 77713.77, 775.52, 1577.34, 7.11, 5573.14, 796.3, 5507.55, 
                 242.73, 529.69, 398.15, 553.57, 5223.69, 0.98),
         SVector(6283.076, 0, 12566.15, 155.42, 3.52, 18849.23, 242.73),
         SVector(0, 6283.08, 12566.15),
         SVector(0.0),
         SVector(84334.662, 5507.553, 5223.69, 2352.87, 1577.34),
         SVector(5507.55, 5223.69),
         SVector(0, 6283.07585, 12566.1517, 77713.7715, 5753.3849, 7860.4194,
                 11506.77, 3930.21, 5884.927, 5507.553, 5223.694, 5573.143, 
                 11790.629, 1577.344, 10977.079, 18849.228, 5486.778, 6069.78,
                 15720.84, 161000.69, 17260.15, 529.69, 83996.85, 71430.7, 
                 2544.31, 775.52, 9437.76, 6275.96, 4694, 8827.39, 19651.05,
                 12139.55, 12036.46, 2942.46, 7084.9, 5088.63, 398.15, 6286.6,
                 6279.55, 10447.39),
         SVector(6283.07585, 12566.1517, 0, 18849.23, 5507.55, 5223.69, 1577.34,
                 10977.08, 6275.96, 5486.78),
         SVector(6283.0758, 12566.152, 0, 77713.77, 5573.14, 18849.23),
         SVector(6283.076, 12566.15),
         SVector(6283.08))

    # Eq. 9 - 10
    @inbounds @inline function calc_L(i, JME)
        sum(A[i][j]*cos(B[i][j] + C[i][j]*JME) for j in 1:length(A[i]))
    end

    # Earth heliocentric longitude (Eq. 11 and 12)
    global heliocentric_longitude
    @inbounds function heliocentric_longitude(JME)
        Ls = SVector{6,Float64}(calc_L(i, JME) for i in 1:6)
        L = (Ls[1] + Ls[2]*JME + Ls[3]*JME^2 + Ls[4]*JME^3 + Ls[5]*JME^4 + 
               Ls[6]*JME^5)/1e8*rad
        L = constrain_angle_circle(L)
        return L
    end

    # Earth heliocentric latitude (paragraph 3.2.7)
    global heliocentric_latitude
    @inbounds function heliocentric_latitude(JME)
        Ls = SVector{2,Float64}(calc_L(i, JME) for i in 7:8)
        L = (Ls[1] + Ls[2]*JME)/1e8*rad
        return L
    end

    # Earth in Astronomical Units (paragraph 3.2.8)
    global heliocentric_radius
    @inbounds function heliocentric_radius(JME)
        Ls = SVector{5,Float64}(calc_L(i, JME) for i in 9:13)
        R = (Ls[1] + Ls[2]*JME + Ls[3]*JME^2 + Ls[4]*JME^3 + Ls[5]*JME^4)/1e8
        return R
    end

end

# Calculate geocentric longitude (Eq. 13)
function geocentric_longitude(JME)
    Θ = heliocentric_longitude(JME) + π*rad
    Θ >= 2π*rad && (Θ = Θ - 2π*rad)
    return Θ
end

# Calculate geocentric latitude (Eq. 14)
function geocentric_latitude(JME)
   β = -heliocentric_latitude(JME)
   return β
end

#=
    Effects of nutation on coordinates
=#

let
    # Table A4.3 on Periodic Terms for the Nutation in Longitude and Obliquity
    Y = (SVector(0, -2, 0, 0, 0, 0, -2, 0, 0, -2, -2, -2, 0, 2, 0, 2, 0, 0, -2, 
                 0, 2, 0, 0, -2, 0, -2, 0, 0, 2, -2, 0, -2, 0, 0, 2, 2, 0, -2, 
                 0, 2, 2, -2, -2, 2, 2, 0, -2, -2, 0, -2, -2, 0, -1, -2, 1, 0, 
                 0, -1, 0, 0, 2, 0, 2),
         SVector(0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
         0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 
         0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -1, 1, -1, -1, 0, -1),
         SVector(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, -1, 0, 1, -1, -1, 1, 2, -2, 
         0, 2, 2, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 2, -1, 1, 0, 1, 0, 0, 1, 2, 1, 
         -2, 0, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0, 1, -2, 1, 1, 1, -1, 3, 0),
         SVector(0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 2, 0, 2, 2, 
                 2, 0, 2, 2, 2, 2, 0, 0, 2, 0, 0, 0, -2, 2, 2, 2, 0, 2, 2, 0, 2, 
                 2, 0, 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 2, 2, 0, 0, 2, 2, 2, 2),
         SVector(1, 2, 2, 2, 0, 0, 2, 1, 2, 2, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1, 2, 
         2, 0, 2, 0, 0, 1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 2, 0, 2, 1, 0, 2, 1, 1, 1, 
         0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 2, 2))

    a = SVector(-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, 
                -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 
                26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, 
                -7, 6, 6, 6, -6, -6, 5, -5, -5, -5, 4, 4, 4, -4, -4, -4, 3, -3, 
                -3, -3, -3, -3, -3, -3)
    b = SVector(-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1, 
                0, 0, 0.1, 0, -0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.1, 0, 0.1, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    c = SVector(92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, 0, -70, -53, 
                0, -33, 26, 32, 27, 0, -24, 16, 13, 0, -12, 0, 0, -10, 0, -8, 7,
                9, 7, 6, 0, 5, 3, -3, 0, 3, 3, 0, -3, -3, 3, 3, 0, 3, 3, 3, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    d = SVector(8.9, -3.1, -0.5, 0.5, -0.1, 0, -0.6, 0, -0.1, 0.3, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0)

    # Mean elongation of the moon from the Sun (Eq. 15)
    calc_xo(JCE) = (297.085036 +  445267.11148JCE -  0.0019142JCE^2 + JCE^3/189474)°
    # Mean anomaly of the Sun (Eq. 16)
    calc_x1(JCE) = (357.52772  + 35999.050340JCE - 0.0001603JCE^2  - JCE^3/300000)°
    # Mean anomaly of the Moon (Eq. 17)
    calc_x2(JCE) = (134.96298 + 477198.867398JCE + 0.0086972JCE^2  + JCE^3/56250)°
    # Moon's argument of latitude (Eq. 18)
    calc_x3(JCE) = (93.27191 + 483202.017538JCE - 0.0036825JCE^2 + JCE^3/327270)°
    # Longitude of the ascending node of the moon's mean orbit (Eq. 19)
    calc_x4(JCE) = (125.04452 - 1934.136261JCE + 0.0020708JCE^2 + JCE^3/450000)°

    # Nutation in longitude (ΔΨ) and ecliptic (Δϵ)
    global nutation
    @inbounds function nutation(JCE)
        X₀ = calc_xo(JCE)
        X₁ = calc_x1(JCE)
        X₂ = calc_x2(JCE)
        X₃ = calc_x3(JCE)
        X₄ = calc_x4(JCE)
        # Eq. 20 & 22
        ΔΨ = sum(((a[j] + b[j]*JCE)*sin(Y[1][j]*X₀ + Y[2][j]*X₁ + Y[3][j]*X₂ + 
                 Y[4][j]*X₃ + Y[5][j]*X₄)/36000000.0)  for j in 1:length(a))°
        # Eq. 21 & 23
        Δϵ = sum(((c[j] + d[j]*JCE)*cos(Y[1][j]*X₀ + Y[2][j]*X₁ + Y[3][j]*X₂ + 
                Y[4][j]*X₃ + Y[5][j]*X₄)/36000000.0) for j in 1:length(c))°
        return (ΔΨ, Δϵ)
    end


end

# Compute true obliquity of the ecliptic
function obliquity(JME, Δϵ)
    # Eq. 24
    U = JME/10
    ϵ₀ = 84381.448 - 4680.93U - 1.55U^2 + 1999.25U^3 - 51.38U^4 - 249.67U^5 -
         39.05U^6 + 7.12U^7 + 27.87U^8 + 5.79U^9 + 2.45U^10
    # Eq. 25
    ϵ = (ϵ₀/3600)° + Δϵ
    return ϵ
end

# Calculate the aberration correction (Eq. 26)
function aberration(R)
    Δτ = (-20.4898/(3600*R))°
    return Δτ
end

# Calculate the apparent sun longitude (Eq. Eq. 27)
function sun_longitude(Θ, ΔΨ, Δτ)
    λ = Θ + ΔΨ + Δτ
    return λ
end

#=
    Apparent sidereal time at Greenwich as a function of julia day and century,
the nutation effect on lontigude and the obliquity of the eliptic
=#
function sidereal_time(JD, JC, ΔΨ, ϵ)
    # Calculate the mean sidereal time at Greenwich (Eq. 28)
    ν₀ = (280.46061837 + 360.98564736629*(JD - 2451545) +
          0.000387933JC^2 - JC^3/38710000)°
    F = ν₀/(360°)
    F = F - trunc(F)
    ν₀ = F < 0.0 ? (1 + F)*360° : F*360°
    # Eq. 29
    ν = ν₀ + ΔΨ*cos(ϵ)
    return ν
end

#=
    Geocentric sun right ascension (Eq. 30)
=#
function geogentric_sun_ascension(λ, ϵ, β)
    α = atan(sin(λ)cos(ϵ) - tan(β)sin(ϵ), cos(λ))
    α = constrain_angle_circle(α)
    return α
end

#=
    Geocentric sun declination (Eq. 31)
=#
function geocentric_sun_declination(λ, ϵ, β)
    δ = asin(sin(β)cos(ϵ) + cos(β)sin(ϵ)sin(λ))
    return δ
end

#=
    Section 3.11 - Observer local hour angle (Eq. 32) as a function of 
    geographical longitude (σ)
=#
function local_hour_angle(σ, ν, α)
    H = ν + σ - α
    H = constrain_angle_circle(H)
    return H
end

#=
    Sections 3.12 - 3.13: Topocentric sun right ascension (α′), sun declination 
(δ′) and local hour angle (H′) as a function of geographical latitude (ϕ) and 
observer elevation (in meters)
=#
function topocentric_quantities(ϕ, E, R, H, δ, α)
    # Equatorial horizontal parallax (Eq. 33)
    ξ = 8.794°/(3600R)
    # Intermediate terms u, x and y (Eqs. 34 - 36)
    u = atan(0.99664719*tan(ϕ))
    x = cos(u) + E/6378140m*cos(ϕ)
    y = 0.99664719sin(u) + E/6378140m*sin(ϕ)
    # Parallax in the sun right ascension (Eq. 37)
    Δα = atan(-x*sin(ξ)*sin(H), cos(δ) - x*sin(ξ)*cos(H))
    # Topocentric right ascension (Eq. 38)
    α′ = α + Δα
    # Topocentric sun declination (Eq. 39)
    δ′ = atan((sin(δ) - y*sin(ξ))*cos(Δα), cos(δ) - x*sin(ξ)cos(H))
    # Topocentric local hour angle
    H′ = H - Δα
    # Return all quantities
    (α′, δ′, H′)
end

#=
    Section 3.14 - Topocentric zenith angle (θ) including effect of refraction
=#
function topocentric_zenith(ϕ, δ′, H′, P, T, R, Rₛᵣ)
    # Topocentric elevation angle in the absence of refration (Eq. 41)
    e₀ = asin(sin(ϕ)*sin(δ′) + cos(ϕ)*cos(δ′)*cos(H′))
    # Refraction correction when sun is above horizon only (Eq. 42)
    if e₀ >= -(0.26667° + Rₛᵣ)
        Δe = (P/1010mbar)*(283K/T)*1.02°/(60tan(e₀ + 10.3°^2/(e₀ + 5.11°)))
    else
        Δe = 0.0°
    end
    # Topocentric elevation angle (Eq. 43)
    e = e₀ + Δe
    # Topocentric zenith angle (Eq. 44)
    θ = (π/2)rad - e
    return θ
end

#=
    Section 3.15 - Topocentric azimuth angles (Φ and Γ)
=#
function topocentric_azimuth(ϕ, H′, δ′)
    # Topocentric astronomers azimuth angle (Eq. 45)
    Γ = atan(sin(H′), cos(H′)sin(ϕ) - tan(δ′)*cos(ϕ))
    Φ = Γ + π*rad
    Γ = constrain_angle_circle(Γ)
    # Topocentric azimuth angle (Eq. 46)
    Φ = constrain_angle_circle(Φ)
    (Γ, Φ)
end


#=
    Section 3.16: Incidence angle on an oriented surface (Eq. 47 parameterized
by the azimuth angle for solar radiation users)
=#
function incidence_angle(θ, Γ, ω, γ)
    acos(cos(θ)cos(ω) + sin(ω)*sin(θ)*cos(Γ - γ))
end


#=
    Section A1: Equation of Time or difference between solar apparent and mean
solar time.
=#
function equation_of_time(JME, α, ΔΨ, ϵ)
    # Eq. A2
    M = (280.4664567 + 360007.6982779JME + 0.03032028JME^2 +
         JME^3/49931 - JME^4/15300 - JME^5/2000000)°
    M = constrain_angle_circle(M)
    # Eq. A1
    Eang = convert(typeof(1.0°), M - 0.0057183° - α + ΔΨ*cos(ϵ))
    E = Eang*4minute/°
    E < -20.0minute && (E =  E + 1440minute)
    E >  20.0minute && (E =  E - 1440minute)
    return E
end

#=
    Section A2: Sunrise, Sun Transit and Sunset.
    Notice that dom should be calculated for 0 UT
=#
function constrain_0_1(m)
    L = m - floor(m)
    L < zero(typeof(L)) && (L = L + one(typeof(m)))
    return L
end
function sunrise_transit_sunset(year, month, dom, ϕ, σ, ΔT, Rₛᵣ, tz)

    h₀′ = -(0.26667° + Rₛᵣ)

    # Convert UT1 in year-month-dom to Julian
    JD = julian_day(year, month, dom)
    JC = JD_to_JC(JD)

    # Geocentric right ascension and declination at different days (Section A.2.2)
    α₀, δ₀, ΔΨ₀, ϵ₀, R = geocentric_coordinates(year, month, dom, 0.0s)
    αₘ, δₘ, ΔΨₘ, ϵₘ, R = geocentric_coordinates(year, month, dom - 1, 0.0s)
    αₚ, δₚ, ΔΨₚ, ϵₚ, R  = geocentric_coordinates(year, month, dom + 1, 0.0s)

    # Apparent sidereal time at Greenwich (Section A.2.1)
    ν = sidereal_time(JD, JC, ΔΨ₀, ϵ₀)

    # Approximate sun transit time (Section A.2.3)
    m₀ = (α₀ - σ - ν)/(2π*rad)
    
    # Local hour angle when sun elevation equals 0.8333° (Eq. A4)
    cosH₀ = (sin(h₀′) - sin(ϕ)sin(δ₀))/(cos(ϕ)cos(δ₀))
    if abs(cosH₀) > 1.0
        return (-99999.0, -99999.0, -99999.0)
    else
        H₀ = constrain_angle_half_circle(acos(cosH₀))
        H₀ < 0.0 && (return (-99999.0, -99999.0, -99999.0))
    end
  
    # Approximate sunrise time as fraction of day (Eq. A5)
    m₁ = m₀ - H₀/(2π*rad)

    # Approximate sunset time as fraction of day (Eq. A6)
    m₂ = m₀ + H₀/(2π*rad)

    # Limit approximate transit, sunrise and sunset times (Section A.2.7)
    m₀ = constrain_0_1(m₀)
    m₁ = constrain_0_1(m₁)
    m₂ = constrain_0_1(m₂)

    # Convert the times above to sidereal time at Greenwich using Earth 
    # sidereal rotation rate (Eq. A7)
    ν₀ = ν + 360.985647°*m₀
    ν₁ = ν + 360.985647°*m₁
    ν₂ = ν + 360.985647°*m₂

    # Convert to terrestrial time (Eq. A8)
    n₀ = JD_to_JDE(m₀, ΔT)
    n₁ = JD_to_JDE(m₁, ΔT)
    n₂ = JD_to_JDE(m₂, ΔT)

    # Calculate α′ᵢ for each time (Eq. A9)
    a = α₀ - αₘ
    a > 2° && (a = constrain_0_1(a))
    b = αₚ - α₀
    b > 2° && (b = constrain_0_1(b))
    c = b - a
    α′₀ = α₀ + n₀*(a + b + c*n₀)/2
    α′₁ = α₀ + n₁*(a + b + c*n₁)/2
    α′₂ = α₀ + n₂*(a + b + c*n₂)/2

    # Calculate δ′ᵢ for each time (Eq. A10)
    a′ = δ₀ - δₘ
    a′ > 2° && (a′ = constrain_0_1(a′))
    b′ = δₚ - δ₀
    b′ > 2° && (b′ = constrain_0_1(b′))
    c′ = b′ - a′
    δ′₀ = δ₀ + n₀*(a′ + b′ + c′*n₀)/2
    δ′₁ = δ₀ + n₁*(a′ + b′ + c′*n₁)/2
    δ′₂ = δ₀ + n₂*(a′ + b′ + c′*n₂)/2

    # Calculate local hour angle for each time (Eq. A11)
    H′₀ = constrain_angle_half_circle_pm(ν₀ + σ - α′₀)
    H′₁ = constrain_angle_half_circle_pm(ν₁ + σ - α′₁)
    H′₂ = constrain_angle_half_circle_pm(ν₂ + σ - α′₂)

    # Calculate sun altitude for each time (Eq. A12)
    h₀ = asin(sin(ϕ)sin(δ′₀) + cos(ϕ)cos(δ′₀)cos(H′₀))
    h₁ = asin(sin(ϕ)sin(δ′₁) + cos(ϕ)cos(δ′₁)cos(H′₁))
    h₂ = asin(sin(ϕ)sin(δ′₂) + cos(ϕ)cos(δ′₂)cos(H′₂))

    # Sun transit as fraction of day (A13) in the UT1 scale
    T = m₀ - H′₀/(2π*rad)

    # Sunrise as fraction of day (A14) in the UT1 scale
    R = m₁ + (h₁ - h₀′)/(2π*rad*cos(δ′₁)cos(ϕ)sin(H′₁))

    # Sunset as fraction of day (A14 as modified by A.2.15) in the UT1 scale
    S = m₂ + (h₂ - h₀′)/(2π*rad*cos(δ′₂)cos(ϕ)sin(H′₂))

    # Calculate timepoints as time of the day in local time
    R = 24*constrain_0_1(R + tz/24)
    T = 24*constrain_0_1(T + tz/24)
    S = 24*constrain_0_1(S + tz/24)


    return (R, T, S)
end

# Function to compute α and δ and also return ΔΨ and ϵ for posterior calculations
# This code will then be reused to compute transit time, sunrise and sunset
function geocentric_coordinates(year, month, dom, ΔT)
    # Convert UTC in year-month-dom to julian and ephimeral time
    JD = julian_day(year, month, dom)
    JC = JD_to_JC(JD)
    # ΔT = UT1_to_TT(year)
    JDE = JD_to_JDE(JD, ΔT)
    JCE = JDE_to_JCE(JDE)
    JME = JCE_to_JME(JCE)

    # Geocentric coordinates and radius
    R = heliocentric_radius(JME)
    Θ = geocentric_longitude(JME)
    β = geocentric_latitude(JME)

    # Corrected geocentric coordinates
    ΔΨ, Δϵ = nutation(JCE)
    Δτ = aberration(R)
    ϵ = obliquity(JME, Δϵ)
    λ = sun_longitude(Θ, ΔΨ, Δτ)

    # Equatorial coordinate system
    α = geogentric_sun_ascension(λ, ϵ, β)
    δ = geocentric_sun_declination(λ, ϵ, β)

    # Return the coordinate system + terms required for sidereal time
    return (α, δ, ΔΨ, ϵ, R)
end

#=
    This function calculates the same output as the official SPA implementation
in the C language. It can be used as a drop-in substitute and for testing this
implementation against the original.

Inputs:
    - Year
    - Month
    - Day 
    - Hour 
    - Minute 
    - Second 
    - ΔUT1: Difference between UTC and UT
    - ΔT: Difference between Earth rotation time and terrestrial dynamical time
    - Timezone: Difference between local time and UTC (negative West of Greenwich)
    - ϕ: Geographical latitude (negative South of Equator).
    - σ: Geographical longitude (negative west of Greenwich)
    - E: Elevation.
    - P: Annual average local pressure.
    - T: Annual average local temperature.
    - γ: Surface azimuth rotation (eastward from North)
    - ω: Surface elevation angle with respect to horizontal plane
    - Rₛᵣ: Atmospheric refraction at sunrise and sunset (default is 0.5667°)

Outputs:
    - θ: Topocentric zenith angle.
    - Φₐ: Topocentric azimuth angle (westward from South).
    - Φ: Topocentric azimuth angle (eastward from North).
    - I: Surface incidence angle.
    - ST: Local sun transit time (i.e. solar noon) as fraction of day.
    - SR: Local sunrise time as fraction of day.
    - SS: Local sunset time as fraction of day.
=#
function SPA(year, month, day, hour, minute, second, timezone, 
             ϕ, σ, ω, γ, E, P, T, ΔT, Rₛᵣ = 0.5667°)

    # Calculate day of month in UTc (assume UT1 = UTC)
    dom = day + hour/24 + minute/60/24 + second/3600/24 - timezone/24

    # Equatorial coordinate system (declination and ascension)
    α, δ, ΔΨ, ϵ, R = geocentric_coordinates(year, month, dom, ΔT)

    # Time of the day in angle
    JD = julian_day(year, month, dom)
    JC = JD_to_JC(JD)
    ν = sidereal_time(JD, JC, ΔΨ, ϵ)
    H = local_hour_angle(σ, ν, α)

    # Topocentric variables
    α′, δ′, H′ = topocentric_quantities(ϕ, E, R, H, δ, α)

    # Calculate the local solar angles including effect of refraction
    θ = topocentric_zenith(ϕ, δ′, H′, P, T, R, Rₛᵣ)
    Γ, Φ = topocentric_azimuth(ϕ, H′, δ′)

    # Calculate angle of incidence given a surface
    I = incidence_angle(θ, Γ, ω, γ)

    # Calculate sun transit, sunrise and sunset
    R, T, S = sunrise_transit_sunset(year, month, day, ϕ, σ, ΔT, Rₛᵣ, timezone)

    # Equation of time
    E = equation_of_time(JCE_to_JME(JDE_to_JCE(JD_to_JDE(JD, ΔT))), α, ΔΨ, ϵ)

    # Return all outputs of SPA
    (θ = θ, Γ = Γ, Φ = Φ, I = I, T = T, R = R, S = S, E = E)

end