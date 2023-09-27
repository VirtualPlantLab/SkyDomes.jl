using Sky
using Test
using StaticArrays

let
    threes = @SVector [1.0, 1.0, 1.0]

    # Discretization of the sky using two methods
    angles = equal_angle_intervals(9, 12)
    solid = equal_solid_angles(9, 12)
    @test solid[end][2] ≈ π / 2
    @test solid[end][4] ≈ 2π

    # Uniform distribution - 1 waveband
    uniform = UniformSky()
    solid_uniform = radiosity(uniform, solid)
    @test solid_uniform[end][2] ≈ π / 2
    @test solid_uniform[end][4] ≈ 2π
    @test all(sum(solid_uniform.I) .≈ [1.0])
    angles_uniform = radiosity(uniform, angles)
    @test all(sum(angles_uniform.I) .≈ [1.0])

    # Uniform distribution - 3 wavebands
    uniform = UniformSky()
    solid_uniform = radiosity(uniform, solid, threes)
    @test solid_uniform[end][2] ≈ π / 2
    @test solid_uniform[end][4] ≈ 2π
    @test all(sum(solid_uniform.I) .≈ threes)
    angles_uniform = radiosity(uniform, angles, threes)
    @test all(sum(angles_uniform.I) .≈ threes)

    # Standard overcast distribution - 1 waveband
    standard = StandardSky()
    solid_standard = radiosity(standard, solid)
    @test all(sum(solid_standard.I) .≈ [1.0])
    angles_standard = radiosity(standard, angles)
    @test all(sum(angles_standard.I) .≈ [1.0])

    # Standard overcast distribution - 3 wavebands
    standard = StandardSky()
    solid_standard = radiosity(standard, solid, threes)
    @test all(sum(solid_standard.I) .≈ [1.0, 1.0, 1.0])
    angles_standard = radiosity(standard, angles, threes)
    @test all(sum(angles_standard.I) .≈ [1.0, 1.0, 1.0])

    # CIE distributions - 1 waveband
    for i in 1:15
        cie = CIE(type = i, θₛ = π / 4, Φₛ = π / 2)
        solid_cie = radiosity(cie, solid)
        @test all(sum(solid_cie.I) .≈ [1.0])
        angles_cie = radiosity(cie, angles)
        @test all(sum(angles_cie.I) .≈ [1.0])
    end

    # CIE distributions - 3 wavebands
    for i in 1:15
        cie = CIE(type = i, θₛ = π / 4, Φₛ = π / 2)
        solid_cie = radiosity(cie, solid, threes)
        @test all(sum(solid_cie.I) .≈ [1.0, 1.0, 1.0])
        angles_cie = radiosity(cie, angles, threes)
        @test all(sum(angles_cie.I) .≈ [1.0, 1.0, 1.0])
    end
end
