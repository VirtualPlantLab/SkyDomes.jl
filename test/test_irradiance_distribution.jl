using Sky
using Test

let
    
    # Discretization of the sky using two methods
    angles = equal_angle_intervals(9, 12);
    solid = equal_solid_angles(9,12);
    @test solid[end][2] ≈ π/2
    @test solid[end][4] ≈ 2π

    # Uniform distribution
    uniform = UniformSky()
    solid_uniform = radiosity(uniform, solid, true);
    @test solid_uniform[end][2] ≈ π/2
    @test solid_uniform[end][4] ≈ 2π
    @test sum(solid_uniform.I) ≈ 1.0
    angles_uniform = radiosity(uniform, angles, true);
    @test sum(angles_uniform.I) ≈ 1.0

    # Standard overcast distribution
    standard = StandardSky()
    solid_standard = radiosity(standard, solid, true);
    @test sum(solid_standard.I) ≈ 1.0
    angles_standard = radiosity(standard, angles, true);
    @test sum(angles_standard.I) ≈ 1.0

    # CIE distributions
    for i in 1:15
        cie = CIE(i, θₛ = π/4, Φₛ = π/2)
        solid_cie = radiosity(cie, solid, true);
        @test sum(solid_cie.I) ≈ 1.0
        angles_cie = radiosity(cie, angles, true);
        @test sum(angles_cie.I) ≈ 1.0
    end

end