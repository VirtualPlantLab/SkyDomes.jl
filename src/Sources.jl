### This file contains public API ###
# equal_angle_intervals
# equal_solid_angles
# radiance
# radiosity
# CIE
# sky
# StandardSky
# UniformSky

#=
    Angular distribution of diffuse radiation in the sky integrated over
discrete sectors.
=#

abstract type RadianceDistribution end

################################################################################
############################ SkySectors Object #################################
################################################################################

#=
    Object containing the output of every discretization method and convenience
methods to work with it as a collection of values.
=#
struct SkySectors{T}
    θₗ::Vector{T}
    θᵤ::Vector{T}
    Φₗ::Vector{T}
    Φᵤ::Vector{T}
end

function getindex(s::SkySectors, i::Int)
    (θₗ = s.θₗ[i], θᵤ = s.θᵤ[i], Φₗ = s.Φₗ[i], Φᵤ = s.Φᵤ[i])
end
length(s::SkySectors) = length(s.θₗ)
lastindex(s::SkySectors) = length(s)
iterate(s::SkySectors, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

################################################################################
########################### Discretization methods #############################
################################################################################

"""
    equal_angle_intervals(ntheta, nphi)

Discretize the sky into `ntheta` zenith rings of `nphi` sectors each assuming
the same angle intervals for each sector (Δθ = 90/ntheta and ΔΦ = 360/nphi).
Returns an object of type `SkySectors`. See package documentation for details.
"""
function equal_angle_intervals(ntheta, nphi)
    Δθ = 90.0 / ntheta
    ΔΦ = 360.0 / nphi
    θₗ = repeat(0.0:Δθ:(90.0 - Δθ), inner = nphi)
    θᵤ = repeat(Δθ:Δθ:90.0, inner = nphi)
    Φₗ = repeat(0.0:ΔΦ:(360.0 - ΔΦ), outer = ntheta)
    Φᵤ = repeat(ΔΦ:ΔΦ:360.0, outer = ntheta)

    SkySectors(θₗ, θᵤ, Φₗ, Φᵤ)
end

"""
    equal_solid_angles (ntheta, nphi)

Discretize the sky into `ntheta` zenith rings and a number of sectors per ring
that is proportional to `sin(θ)`. The total number of sectors will be `ntheta*nphi`.
Returns an object of type `SkySectors`. See package documentation for details.
"""
function equal_solid_angles(ntheta, nphi)
    # Distribution sectors along zenith angles
    Δθ = 90.0 / ntheta
    uθₗ = 0.0:Δθ:(90.0 - Δθ)
    uθᵤ = Δθ:Δθ:90.0

    # Calculate number of azimuth sectors and ΔΦ per zenith ring
    fac = cosd.(uθₗ) - cosd.(uθᵤ)
    n = ntheta * nphi
    nphis = round.(fac ./ sum(fac) .* n)
    sum(nphis) != n && (nphis[end] = nphis[end] + n - sum(nphis))
    ΔΦs = 360.0 ./ nphis

    # Generate coordinates of all sectors
    c = 1
    θₗ = Vector{typeof(Δθ)}(undef, n)
    θᵤ = Vector{typeof(Δθ)}(undef, n)
    Φₗ = Vector{typeof(Δθ)}(undef, n)
    Φᵤ = Vector{typeof(Δθ)}(undef, n)
    for i in 1:ntheta
        ΔΦ = ΔΦs[i]
        for j in 1:nphis[i]
            θₗ[c] = Δθ * (i - 1)
            θᵤ[c] = Δθ * i
            Φₗ[c] = ΔΦ * (j - 1)
            Φᵤ[c] = ΔΦ * j
            c += 1
        end
    end

    SkySectors(θₗ, θᵤ, Φₗ, Φᵤ)
end

################################################################################
################################## Sky dome ####################################
################################################################################

#=
    Object containing sky sectors and the radiosity associated to each of them
=#
struct SkyDome{AT, IT}
    sectors::SkySectors{AT}
    I::Vector{IT}
end

function getindex(s::SkyDome, i::Int)
    (θₗ = s.sectors.θₗ[i], θᵤ = s.sectors.θᵤ[i],
        Φₗ = s.sectors.Φₗ[i], WΦᵤ = s.sectors.Φᵤ[i],
        I = s.I[i])
end
length(s::SkyDome) = length(s.sectors.θₗ)
lastindex(s::SkyDome) = length(s)
iterate(s::SkyDome, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

function DirectionalSource(sky::SkyDome, box, nrays; α = 180.0, alpha_soil = 0.0, beta_soil = 180.0)
    # Determine the zenith and azimuth angle of each sector
    θₗ = convert(Vector{Float64}, getproperty.(sky.sectors, :θₗ))
    θᵤ = convert(Vector{Float64}, getproperty.(sky.sectors, :θᵤ))
    Φₗ = convert(Vector{Float64}, getproperty.(sky.sectors, :Φₗ))
    Φᵤ = convert(Vector{Float64}, getproperty.(sky.sectors, :Φᵤ))
    θ = (θₗ .+ θᵤ) ./ 2
    Φ = (Φₗ .+ Φᵤ) ./ 2
    # Distribute the total number of rays across the sources
    nraysi = Int.(round.(nrays / length(θₗ)))
    # Create an array of directional sources
    [DirectionalSource(box, θ = θ[i], Φ = Φ[i], α = α, alpha_soil = alpha_soil,
        beta_soil = beta_soil, radiosity = sky.I[i], nrays = nraysi) for i in eachindex(θ)]
end

#This function creates a sky dome of diffuse irradiance for a given mesh, using
# different models of angular distribution (sky_model) and methods of
# discretization (dome_method). It returns a vector of directional sources as
# required by the ray tracer in VPL.
function sky_dome(box; Idif = 1.0, nrays_dif = 1_000,
    sky_model = StandardSky, dome_method = equal_solid_angles, ntheta = 9,
    nphi = 12, α = 180.0, alpha_soil = 0.0, beta_soil = 180.0, kwargs...)
    # Generate the angular distribution of irradiance according to different methods
    if sky_model == CIE
        sky_distro = CIE(kwargs...)
    else
        sky_distro = sky_model()
    end
    # Discretization method
    dome_mesh = dome_method(ntheta, nphi)
    # Generate the sky sectors and their radiosity
    skydome = radiosity(sky_distro, dome_mesh, Idif)
    # Convert the sky dome into a vector of DirectionalSources
    sources = DirectionalSource(skydome, box, nrays_dif; α = α,
        alpha_soil = alpha_soil, beta_soil = beta_soil)
    return sources
end

################################################################################
#################################### Sky  ######################################
################################################################################

function sky(mesh::Mesh; kwargs...)
    @error "Sky dome can only be used with a grid cloner. Please run the `accelerate` function on the mesh and make sure a grid cloner is created."
end

"""
    sky(mesh; Idir = 0.77, nrays_dir = 100_000, theta_dir = 0.0, phi_dir = 0.0,
               Idif = 0.23, nrays_dif = 1_000_000, sky_model = StandardSky,
               dome_method = equal_solid_angles, ntheta = 9, nphi = 12,
               α = 180.0, alpha_soil = 0.0, beta_soil = 180.0, warn_slope = true,
               kwargs...)

Create a vector of directional radiation sources representing diffuse and
direct solar radiation for a given mesh.

The soil surface is assumed to coincide with the XY plane. The `alpha_soil` and
`beta_soil` angles describe a sloped soil and are applied to both the direct and the
diffuse radiation, tilting the incoming directions relative to this fixed soil frame.
Light sources whose rays fall below the horizon of the sloped soil cannot reach the
surface and are removed (a warning is emitted once, see `warn_slope`).

# Arguments
- `mesh`: A `Mesh` object generated by VPL.
- `Idir`: The direct solar radiation measured on the horizontal plane (a single value or tuple).
- `nrays_dir`: The number of rays to be generated for direct solar radiation.
- `theta_dir`: The zenith angle of the sun position (degrees).
- `phi_dir`: The azimuthal angle of the sun position (degrees).
- `Idif`: The diffuse solar radiation measured on the horizontal plane (a single value or tuple).
- `nrays_dif`: The total number of rays to be generated diffuse solar radiation.
- `sky_model`: The angular distribution of diffuse irradiance (`StandardSky`, `UniformSky` or `CIE`).
- `dome_method`: The method to discretize hemisphere into patches for diffuse solar radiation (`equal_solid_angles` or `equal_angle_intervals`).
- `ntheta`: The number of divisions along the zenith angle for `dome_method`.
- `nphi`: The number of divisions along the azimuthal angle for `dome_method`.
- `α`: The azimuth of the local X axis (degrees, default `180.0`).
- `alpha_soil`: The slope inclination of the soil (degrees; `0` = horizontal, `90` = vertical; default `0`).
- `beta_soil`: The azimuth of the slope normal (degrees; e.g. `180` = south-facing; default `180`).
- `warn_slope`: Whether to emit a warning when sources are removed because they fall below a sloped soil's horizon (default `true`).
- `kwargs...`: Additional arguments to be used when `dome_method = CIE`

# Returns
A vector of directional sources that can be used for ray tracing calculations in VPL.
```
"""
function sky(mesh::AccMesh;
    # Inputs for direct solar radiation
    Idir = 0.77,
    nrays_dir = 100_000,
    theta_dir = 0.0,
    phi_dir = 0.0,
    # Inputs for diffuse solar radiation
    Idif = 0.23,
    nrays_dif = 1_000_000,
    sky_model = StandardSky,
    dome_method = equal_solid_angles,
    ntheta = 9,
    nphi = 12,
    # Orientation of the local coordinate system and (possibly sloped) soil
    α = 180.0,
    alpha_soil = 0.0,
    beta_soil = 180.0,
    warn_slope = true,
    kwargs...)
    # Ensure valid inputs
    any(Idir .< 0.0) && @warn "Idir should be non-negative. Proceed at your own risk."
    any(Idif .< 0.0) && @warn "Idif should be non-negative. Proceed at your own risk."
    mesh.grid.nleaves <= 1 && @warn "The scene should include clones in the x and y directions when using Sky domes. Proceed at your own risk."
    # To avoid multiple repetitions
    has_diffuse = false
    has_direct = false
    # Use the correct box to account for overlap in the grid cloner
    box_orig = mesh.acc.gbox
    c        = center(box_orig)
    boxmin   = Vec(c[1] - mesh.grid.dx/2, c[2] - mesh.grid.dy/2, box_orig.min[3])
    boxmax   = Vec(c[1] + mesh.grid.dx/2, c[2] + mesh.grid.dy/2, box_orig.max[3])
    box      = AABB(boxmin, boxmax)
    # Generation dome with diffuse solar radiation
    if any(Idif .> 0.0)
        @assert ntheta > 0
        @assert nphi > 0
        sources = sky_dome(box, Idif = convert_svector(Idif), nrays_dif = nrays_dif,
            sky_model = sky_model, dome_method = dome_method,
            ntheta = ntheta, nphi = nphi, α = α, alpha_soil = alpha_soil,
            beta_soil = beta_soil, kwargs...)
        has_diffuse = true
    end
    # Generate directional source for direct solar radiation
    if any(Idir .> 0.0)
        source = DirectionalSource(box, θ = theta_dir, Φ = phi_dir, α = α,
                                   alpha_soil = alpha_soil, beta_soil = beta_soil,
                                   radiosity = convert_svector(Idir), nrays = nrays_dir)
        has_direct = true
        has_diffuse && push!(sources, source)
    end
    # Assemble the full set of sources
    if has_diffuse
        result = sources
    elseif has_direct
        result = [source]
    else
        error("Attempt to create a sky with Idir = Idif = 0")
    end
    # Remove sources that fall below the soil horizon and cannot reach the surface.
    # The soil is the XY plane (normal = (0,0,1)); the slope angles tilt the incoming
    # directions. A source reaches the soil only when its ray travels downward.
    reaches(s) = s.angle.dir[3] < 0
    kept = filter(reaches, result)
    n_removed = length(result) - length(kept)
    warn_slope && n_removed > 0 && @warn "$(n_removed) light source(s) removed because they fall below the horizon of the sloped soil (alpha_soil = $alpha_soil, beta_soil = $beta_soil) and cannot reach the surface. Pass warn_slope = false to silence."
    isempty(kept) && error("All light sources fall below the soil horizon; no radiation reaches the surface.")
    return kept
end

# Make sure that irradiance is converted to the right format (SVector)
function convert_svector(I)
    if I isa Number
        return SVector{1, Float64}(I)
    elseif I isa Tuple
        return SVector{length(I), Float64}(I)
    else
        return I
    end
end

################################################################################
############################## Uniform radiation ###############################
################################################################################

"""
    UniformSky()

Model of uniform sky diffuse radiation. See package documentation for details.
"""
struct UniformSky <: RadianceDistribution
end

"""
    radiosity(m::UniformSky, sky::SkySectors, Idif::SVector{nw, Float64}) where nw

Calculate the radiosity of each section of `sky` on the horizontal plane given diffuse
irradiance on the horizontal plane (`Idif` with `nw` wavebands) assuming a Uniform Sky model.
See package documentation for details.
"""
function radiosity(m::UniformSky,
    sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0)) where {nw}
    # Equation 15
    I = SVector{nw, Float64}[(t = (s.Φᵤ - s.Φₗ) * (cosd(s.θₗ)^2 - cosd(s.θᵤ)^2) / 2 / 180;
                              SVector{nw, Float64}(t .* Idif[i] for i in 1:nw)) for s in sky]
    SkyDome(sky, I)
end

################################################################################
############################# Standard overcast sky ############################
################################################################################

"""
    StandardSky()

Standard model of overcast sky diffuse radiation Moon & Spencer (1942).
See package documentation for details.
"""
struct StandardSky <: RadianceDistribution
end

"""
    radiosity(m::StandardSky, sky::SkySectors, Idif::SVector{nw, Float64})

Calculate the radiosity of each section of `sky` on the horizontal plane given diffuse
irradiance on the horizontal plane (`Idif` with `nw` wavebands) assuming a Standard Sky
model and for `nw` wavebands. See package documentation for details.
"""
function radiosity(m::StandardSky,
    sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0)) where {nw}
    # Equation 22
    I = SVector{nw, Float64}[(t = (s.Φᵤ - s.Φₗ) * ((4cosd(s.θₗ) + 3) * cosd(s.θₗ)^2 -
                                   (4cosd(s.θᵤ) + 3) * cosd(s.θᵤ)^2) / 14 / 180;
    SVector{nw, Float64}(t .* Idif[i] for i in 1:nw)) for s in sky]
    SkyDome(sky, I)
end

################################################################################
############################ General CIE sky models ############################
################################################################################

#=
    General Sky model from CIE Standard (Darula & Kittler, 2002)
=#
struct CIE <: RadianceDistribution
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    θₛ::Float64
    Φₛ::Float64
    denom::Float64
end

"""
    CIE(;type = 1, θₛ = 0.0, Φₛ = 0.0, rtol = sqrt(eps(Float64)), atol = 0.0,
maxevals = typemax(Int))

Create a standard CIE model of sky diffuse radiance on the horizontal plane as described by
Darula and Kittler (2002). The argument `type` can have values from 1 to 15
representing the 15 standard CIE models. θₛ and Φₛ are the zenith and azimuth
angles of the solar disc in degrees. `rtol` and `atol` and `maxevals` are the relative
tolerance, absolute tolerance and maximum number of function evaluation of the
numerical integration algorithm. See package documentation for details.
"""
function CIE(; type = 1, θₛ = 0.0, Φₛ = 0.0, rtol = sqrt(eps(Float64)), atol = 0.0,
    maxevals = typemax(Int))
    a, b, c, d, e = CIEmodel(type)
    denom = hcubature(x -> f_denom(x[1], x[2], a, b, c, d, e, θₛ, Φₛ),
        (0.0, 0.0), (90.0, 360.0); rtol = rtol, atol = atol,
        maxevals = maxevals)
    CIE(a, b, c, d, e, θₛ, Φₛ, denom[1])
end

# Generate the parameters for the 15 CIE models
function CIEmodel(type::Int)
    if type == 1
        (a = 4.0, b = -0.7, c = 0.0, d = -1.0, e = 0.0)
    elseif type == 2
        (a = 4.0, b = -0.7, c = 2.0, d = -1.5, e = 0.15)
    elseif type == 3
        (a = 1.1, b = -0.8, c = 0.0, d = -1.0, e = 0.0)
    elseif type == 4
        (a = 1.1, b = -0.8, c = 2.0, d = -1.5, e = 0.15)
    elseif type == 5
        (a = 0.0, b = -1.0, c = 0.0, d = -1.0, e = 0.0)
    elseif type == 6
        (a = 0.0, b = -1.0, c = 2.0, d = -1.5, e = 0.15)
    elseif type == 7
        (a = 0.0, b = -1.0, c = 5.0, d = -2.5, e = 0.30)
    elseif type == 8
        (a = 0.0, b = -1.0, c = 10.0, d = -3.0, e = 0.45)
    elseif type == 9
        (a = -1.0, b = -0.55, c = 2.0, d = -1.5, e = 0.15)
    elseif type == 10
        (a = -1.0, b = -0.55, c = 5.0, d = -2.5, e = 0.30)
    elseif type == 11
        (a = -1.0, b = -0.55, c = 10.0, d = -3.0, e = 0.45)
    elseif type == 12
        (a = -1.0, b = -0.32, c = 10.0, d = -3.0, e = 0.45)
    elseif type == 13
        (a = -1.0, b = -0.32, c = 16.0, d = -3.0, e = 0.45)
    elseif type == 14
        (a = -1.0, b = -0.15, c = 16.0, d = -3.0, e = 0.3)
    elseif type == 15
        (a = -1.0, b = -0.15, c = 24.0, d = -2.8, e = 0.15)
    else
        error("The type of CIE model must be indicated with a number between 1 and 15")
    end
end

# Function used to compute the denom normalisation factor (θ, Φ, θₛ, Φₛ in degrees)
function f_denom(θ, Φ, a, b, c, d, e, θₛ, Φₛ)
    # γₛ is kept in radians: it appears in exp(d*γₛ) where d was calibrated for radian angles
    # Note: tiny numerical errors introduced in trigonometric function, hence clamps
    γₛ = acos(clamp(cosd(θₛ) * cosd(θ) + sind(θₛ) * sind(θ) * cosd(abs(Φ - Φₛ)), -1.0, 1.0))
    fCIE(a, b, c, d, e, θ, γₛ) * sind(θ) * cosd(θ)
end

# Relative radiance with respect to any direction in the sky dome (θ in degrees, γₛ in radians)
function fCIE(a, b, c, d, e, θ, γₛ)
    # γₛ stays in radians here: exp(d*γₛ) and exp(d*π/2) use the CIE parameters as calibrated
    (1 + c * (exp(d * γₛ) - exp(d * π / 2)) + e * cos(γₛ)^2) * (1 + a * exp(b / cosd(θ)))
end

# Radiance for a given angle normalized by horizontal irradiance (θ, Φ in degrees)
function radiance(m::CIE, θ, Φ)
    # γₛ is kept in radians: it appears in exp(d*γₛ) where d was calibrated for radian angles
    # Note: tiny numerical errors introduced in trigonometric function, hence clamps
    γₛ = acos(clamp(cosd(m.θₛ) * cosd(θ) + sind(m.θₛ) * sind(θ) * cosd(abs(Φ - m.Φₛ)), -1.0, 1.0))
    fCIE(m.a, m.b, m.c, m.d, m.e, θ, γₛ) / m.denom
end

# Radiosity for a sector of the sky on horizontal plane normalized by horizontal irradiance
function radiosity(m::CIE, sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0);
    rtol = sqrt(eps(Float64)), atol = 0.0, maxevals = typemax(Int)) where {nw}
    I = SVector{nw, Float64}[(t = hcubature(x -> radiance(m, x[1], x[2]) * cosd(x[1]) *
                                                 sind(x[1]),
        (s.θₗ, s.Φₗ), (s.θᵤ, s.Φᵤ); rtol = rtol, atol = atol,
        maxevals = maxevals)[1];
    SVector{nw, Float64}(t .* Idif[i] for i in 1:nw)) for s in sky]
    SkyDome(sky, I)
end
