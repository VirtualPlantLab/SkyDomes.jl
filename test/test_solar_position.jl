using Test
import Sky
S = Sky
import Unitful: rad

let 

# Test calculation of julian day from universal time
@test S.julian_day(2000, 1, 1.5) == 2451545.0
@test S.julian_day(1999, 1, 1.0) == 2451179.5
@test S.julian_day(1987, 1, 27) == 2446822.5
@test S.julian_day(1987, 6, 19.5) == 2446966.0
@test S.julian_day(1988, 1, 27) == 2447187.5
@test S.julian_day(1988, 6, 19.5) == 2447332.0
@test S.julian_day(1900, 1, 1.0)  == 2415020.5
@test S.julian_day(1600, 1, 1.0) == 2305447.5
@test S.julian_day(1600, 12, 31) == 2305812.5
@test S.julian_day(837, 4, 10 + 7/24 + 12/60/24, false) == 2026871.8
@test S.julian_day(-123, 12, 31, false) == 1676496.5
@test S.julian_day(-122, 1, 1, false) == 1676497.5
@test S.julian_day(-1000, 7, 12.5, false) == 1356001.0
@test S.julian_day(-1001, 8, 17 + 21/24 + 36/60/24, false) == 1355671.4
@test S.julian_day(-4712, 1, 1.5, false) == 0.0

# Test the correction term between universal time and terrestrial dynamical time
# <<-- These tests are not longer valid as the code is being refactored.... -->>
# @test abs(Sky.UTC_to_TDT(-500) - 17190) < 430
# @test abs(Sky.UTC_to_TDT(-400) - 15530) < 390
# @test abs(Sky.UTC_to_TDT(-300) - 14080) < 360
# @test abs(Sky.UTC_to_TDT(-200) - 12790) < 330
# @test abs(Sky.UTC_to_TDT(-100) - 11640) < 290
# @test abs(Sky.UTC_to_TDT(0) - 10580)    < 260
# @test abs(Sky.UTC_to_TDT(100) - 9600)   < 240
# @test abs(Sky.UTC_to_TDT(200) - 8640)   < 210
# @test abs(Sky.UTC_to_TDT(300) - 7680)   < 180
# @test abs(Sky.UTC_to_TDT(400) - 6700)   < 160
# @test abs(Sky.UTC_to_TDT(500) - 5710)   < 140
# @test abs(Sky.UTC_to_TDT(600) - 4740)   < 120
# @test abs(Sky.UTC_to_TDT(700) - 3810)   < 100
# @test abs(Sky.UTC_to_TDT(800) - 2960)   < 80
# @test abs(Sky.UTC_to_TDT(900) - 2200)   < 70
# @test abs(Sky.UTC_to_TDT(1000) - 1570)  < 55
# @test abs(Sky.UTC_to_TDT(1100) - 1090)  < 40
# @test abs(Sky.UTC_to_TDT(1200) - 740)   < 30
# @test abs(Sky.UTC_to_TDT(1300) - 490)   < 20
# @test abs(Sky.UTC_to_TDT(1400) - 320)   < 20
# @test abs(Sky.UTC_to_TDT(1500) - 200)   < 20
# @test abs(Sky.UTC_to_TDT(1600) - 120)   < 20
# @test abs(Sky.UTC_to_TDT(1700) - 9)     < 5
# @test abs(Sky.UTC_to_TDT(1750) - 13)    < 2
# @test abs(Sky.UTC_to_TDT(1800) - 14)    < 1
# @test abs(Sky.UTC_to_TDT(1850) - 7)     < 1
# @test abs(Sky.UTC_to_TDT(1900) - -3)    < 1
# @test abs(Sky.UTC_to_TDT(1950) - 29)    < 0.1

# Test calculation of month and day of month from day of year
all_days = [Sky.month_day(i)[2] for i in 1:365]
all([count(i -> i == j, all_days) for j = 1:12] .==
    [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

all_days = [Sky.month_day(i, true)[2] for i in 1:366]
all([count(i -> i == j, all_days) for j = 1:12] .==
    [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    
# Test the implementation of SPA (as described in documentation)

end
