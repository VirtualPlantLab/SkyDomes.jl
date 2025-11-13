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

"""
    SkySectors

Data structure that contain the zenith and azimuth angles of the lower and upper
bounds of each sector of the sky dome. The angles are stored in radians.
See package documentation for details.

This type is immutable.

## Fields:
- `θₗ`: Vector of lower zenith angles (in radians) for each sector.
- `θᵤ`: Vector of upper zenith angles (in radians) for each sector.
- `Φₗ`: Vector of lower azimuth angles (in radians) for each sector.
- `Φᵤ`: Vector of upper azimuth angles (in radians) for each sector.
"""
struct SkySectors{T}
    θₗ::Vector{T}
    θᵤ::Vector{T}
    Φₗ::Vector{T}
    Φᵤ::Vector{T}
end

"""
    getindex(s::SkySectors, i::Int)

Get the lower and upper zenith and azimuth angles of the `i`-th sector in the SkySectors object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkySectors`.
- `i`: An integer index of the sector.

## Returns
A tuple containing the lower zenith angle (`θₗ`), upper zenith angle (`θᵤ`),
lower azimuth angle (`Φₗ`), and upper azimuth angle (`Φᵤ`) of the `i`-th sector.
"""
function getindex(s::SkySectors, i::Int)
    (θₗ = s.θₗ[i], θᵤ = s.θᵤ[i], Φₗ = s.Φₗ[i], Φᵤ = s.Φᵤ[i])
end

"""
    length(s::SkySectors)

Get the number of sectors in the SkySectors object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkySectors`.

## Returns
The number of sectors in the SkySectors object.
"""
length(s::SkySectors) = length(s.θₗ)

"""
    lastindex(s::SkySectors)

Get the last index of the sectors in the SkySectors object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkySectors`.

## Returns
The last index of the sectors in the SkySectors object.
"""
lastindex(s::SkySectors) = length(s)

"""
    iterate(s::SkySectors, state = 1)

Iterate over the sectors in the SkySectors object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkySectors`.
- `state`: The current state of the iteration (default is 1).

## Returns
A tuple containing the sector data and the next state, or `nothing` if the
end of the sectors is reached.
"""
iterate(s::SkySectors, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

################################################################################
########################### Discretization methods #############################
################################################################################

"""
    equal_angle_intervals(ntheta, nphi)

Discretize the sky into `ntheta` zenith rings of `nphi` sectors each assuming
the same angle intervals for each sector (Δθ = π/2/ntheta and ΔΦ = 2π/nphi).
See package documentation for details.

## Arguments
- `ntheta`: Number of zenith rings.
- `nphi`: Number of sectors per zenith ring.

## Returns
An object of type `SkySectors`. See package documentation for details.
"""
function equal_angle_intervals(ntheta, nphi)
    Δθ = π / 2 / ntheta * 1.0
    ΔΦ = 2π / nphi * 1.0
    θₗ = repeat(0.0:Δθ:((π / 2) - Δθ), inner = nphi)
    θᵤ = repeat(Δθ:Δθ:(π / 2), inner = nphi)
    Φₗ = repeat(0.0:ΔΦ:((2π) - ΔΦ), outer = ntheta)
    Φᵤ = repeat(ΔΦ:ΔΦ:(2π), outer = ntheta)

    SkySectors(θₗ, θᵤ, Φₗ, Φᵤ)
end

"""
    equal_solid_angles (ntheta, nphi)

Discretize the sky into `ntheta` zenith rings and a number of sectors per ring
that is proportional to `sin(θ)`. The total number of sectors will be `ntheta*nphi`.
See package documentation for details.

## Arguments
- `ntheta`: Number of zenith rings.
- `nphi`: Total number of sectors in the sky dome.

## Returns
An object of type `SkySectors`. See package documentation for details.
"""
function equal_solid_angles(ntheta, nphi)
    # Distribution sectors along zenith angles
    Δθ = (π / 2 / ntheta)
    uθₗ = 0.0:Δθ:((π / 2) - Δθ)
    uθᵤ = Δθ:Δθ:(π / 2)

    # Calculate number of azimuth sectors and ΔΦ per zenith ring
    fac = cos.(uθₗ) - cos.(uθᵤ)
    n = ntheta * nphi
    nphis = round.(fac ./ sum(fac) .* n)
    sum(nphis) != n && (nphis[end] = nphis[end] + n - sum(nphis))
    ΔΦs = (2π ./ nphis)

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

"""
    SkyDome(sectors::SkySectors, I::Vector{IT})

Data structure that represents a sky dome with sectors and their associated radiosity.
It contains the sectors of the sky and the radiosity values for each sector.
See package documentation for details.

This type is immutable.

## Fields:
- `sectors`: An instance of `SkySectors` containing the zenith and azimuth angles of each sector.
- `I`: A vector of radiosity values associated with each sector (in W/m²/sr). 
"""
struct SkyDome{AT, IT}
    sectors::SkySectors{AT}
    I::Vector{IT}
end

"""
    getindex(s::SkyDome, i::Int)

Get the lower and upper zenith and azimuth angles of the `i`-th sector in the SkyDome object,
along with the associated radiosity. See package documentation for details.

## Arguments
- `s`: An instance of `SkyDome`.
- `i`: An integer index of the sector.

## Returns
A tuple containing the lower zenith angle (`θₗ`), upper zenith angle (`θᵤ`),
lower azimuth angle (`Φₗ`), upper azimuth angle (`Φᵤ`) in radians, and 
radiosity (`I`) of the `i`-th sector in W/m²/sr.
"""
function getindex(s::SkyDome, i::Int)
    (θₗ = s.sectors.θₗ[i], θᵤ = s.sectors.θᵤ[i],
        Φₗ = s.sectors.Φₗ[i], WΦᵤ = s.sectors.Φᵤ[i],
        I = s.I[i])
end

"""
    length(s::SkyDome)

Get the number of sectors in the SkyDome object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkyDome`.

## Returns
The number of sectors in the SkyDome object.
"""
length(s::SkyDome) = length(s.sectors.θₗ)

"""
    lastindex(s::SkyDome)

Get the last index of the sectors in the SkyDome object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkyDome`.

## Returns
The last index of the sectors in the SkyDome object.
"""
lastindex(s::SkyDome) = length(s)

"""
    iterate(s::SkyDome, state = 1)

Iterate over the sectors in the SkyDome object.
See package documentation for details.

## Arguments
- `s`: An instance of `SkyDome`.
- `state`: The current state of the iteration (default is 1).

## Returns
A tuple containing the sector data and the next state, or `nothing` if the
end of the sectors is reached.
"""
iterate(s::SkyDome, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

"""
    DirectionalSource(sky, box, nrays)

Create a vector of directional sources from a sky dome object.
See package documentation for details.

## Arguments
- `sky`: An instance of `SkyDome` containing the sectors and their radiosity.
- `box`: An axis-aligned bounding box that defines the spatial extent of the sources.
- `nrays`: The total number of rays to be generated for each sector.

## Returns
A vector of `DirectionalSource` objects, each representing a directional source
with its zenith angle (`θ`), azimuth angle (`Φ`), radiosity (`radiosity`), and the 
number of rays (`nraysi`) to be generated.
"""
function DirectionalSource(sky::SkyDome, box, nrays)
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
    [DirectionalSource(box, θ = θ[i], Φ = Φ[i], radiosity = sky.I[i],
        nrays = nraysi) for i in eachindex(θ)]
end

"""
    sky_dome(box; Idif = 1.0, nrays_dif = 1_000,
        sky_model = StandardSky, dome_method = equal_solid_angles,
        ntheta = 9, nphi = 12, kwargs...)
    
Create a sky dome of diffuse irradiance for a given mesh, using
different models of angular distribution (`sky_model`) and methods of
discretization (`dome_method`).

## Arguments
- `box`: An axis-aligned bounding box (AABB) that defines the spatial extent of the sources.
- `Idif`: The diffuse solar radiation measured on the horizontal plane (a single value or tuple).
- `nrays_dif`: The total number of rays to be generated for diffuse solar radiation.
- `sky_model`: The angular distribution of diffuse irradiance (`StandardSky`, `UniformSky` or `CIE`).
- `dome_method`: The method to discretize hemisphere into patches for diffuse solar radiation 
(`equal_solid_angles` or `equal_angle_intervals`).
- `ntheta`: The number of divisions along the zenith angle for `dome_method`.
- `nphi`: The number of divisions along the azimuthal angle for `dome_method`.
- `kwargs...`: Additional arguments to be used when `dome_method = CIE`

## Returns
A vector of `DirectionalSource` objects as required by the ray tracer in VPL.
"""
function sky_dome(box; Idif = 1.0, nrays_dif = 1_000,
    sky_model = StandardSky, dome_method = equal_solid_angles, ntheta = 9,
    nphi = 12, kwargs...)
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
    sources = DirectionalSource(skydome, box, nrays_dif)
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
               kwargs...)

Create a vector of directional radiation sources representing diffuse and
direct solar radiation for a given mesh. See package documentation for details.

# Arguments
- `mesh`: A `Mesh` object generated by VPL.
- `Idir`: The direct solar radiation measured on the horizontal plane (a single value or tuple).
- `nrays_dir`: The number of rays to be generated for direct solar radiation.
- `theta_dir`: The zenith angle of the sun position (radians).
- `phi_dir`: The azimuthal angle of the sun position (radians).
- `Idif`: The diffuse solar radiation measured on the horizontal plane (a single value or tuple).
- `nrays_dif`: The total number of rays to be generated diffuse solar radiation.
- `sky_model`: The angular distribution of diffuse irradiance (`StandardSky`, `UniformSky` or `CIE`).
- `dome_method`: The method to discretize hemisphere into patches for diffuse solar radiation (`equal_solid_angles` or `equal_angle_intervals`).
- `ntheta`: The number of divisions along the zenith angle for `dome_method`.
- `nphi`: The number of divisions along the azimuthal angle for `dome_method`.
- `kwargs...`: Additional arguments to be used when `dome_method = CIE`

# Returns
A vector of directional sources that can be used for ray tracing calculations in VPL.
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
            ntheta = ntheta, nphi = nphi, kwargs...)
        has_diffuse = true
    end
    # Generate directional source for direct solar radiation
    if any(Idir .> 0.0)
        source = DirectionalSource(box, θ = theta_dir, Φ = phi_dir,
                                   radiosity = convert_svector(Idir), nrays = nrays_dir)
        has_direct = true
        has_diffuse && push!(sources, source)
    end
    # Return correct results as a vector of sources
    if has_diffuse
        return sources
    elseif has_direct
        return [source]
    else
        error("Attempt to create a sky with Idir = Idif = 0")
    end
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

Data structure representing a uniform sky diffuse radiation model.
See package documentation for details.

This type is immutable and does not contain any fields.
"""
struct UniformSky <: RadianceDistribution
end

"""
    radiosity(m::UniformSky, sky::SkySectors, Idif::SVector{nw, Float64}) where nw

Calculate the radiosity of each section of `sky` on the horizontal plane given diffuse
irradiance on the horizontal plane (`Idif` with `nw` wavebands) assuming a Uniform Sky model.
See package documentation for details.

## Arguments
- `m`: An instance of `UniformSky`.
- `sky`: An instance of `SkySectors` containing the zenith and azimuth angles of each sector.
- `Idif`: A vector of diffuse irradiance values on the horizontal plane with `nw` wavebands.

## Returns
An instance of `SkyDome` containing the sectors and their associated radiosity values.
"""
function radiosity(m::UniformSky,
    sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0)) where {nw}
    # Equation 15
    I = SVector{nw, Float64}[(t = (s.Φᵤ - s.Φₗ) * (cos(s.θₗ)^2 - cos(s.θᵤ)^2) / 2 / π;
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

## Arguments
- `m`: An instance of `StandardSky`.
- `sky`: An instance of `SkySectors` containing the zenith and azimuth angles of each sector.
- `Idif`: A vector of diffuse irradiance values on the horizontal plane with `nw` wavebands, in W/m²/sr.

## Returns
An instance of `SkyDome` containing the sectors and their associated radiosity values.
"""
function radiosity(m::StandardSky,
    sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0)) where {nw}
    # Equation 22
    I = SVector{nw, Float64}[(t = (s.Φᵤ - s.Φₗ) * ((4cos(s.θₗ) + 3) * cos(s.θₗ)^2 -
                                   (4cos(s.θᵤ) + 3) * cos(s.θᵤ)^2) / 14 / π;
    SVector{nw, Float64}(t .* Idif[i] for i in 1:nw)) for s in sky]
    SkyDome(sky, I)
end

################################################################################
############################ General CIE sky models ############################
################################################################################

#=
    General Sky model from CIE Standard (Darula & Kittler, 2002)
=#
"""
    CIE

Data structure representing a CIE standard sky diffuse radiation model (Darula & Kittler, 2002).
See package documentation for details.

This type is immutable and contains the parameters of the CIE model.

## Fields:
- `a`: Parameter a of the CIE model.
- `b`: Parameter b of the CIE model.
- `c`: Parameter c of the CIE model.
- `d`: Parameter d of the CIE model.
- `e`: Parameter e of the CIE model.
- `θₛ`: Zenith angle of the solar disc (in radians).
- `Φₛ`: Azimuth angle of the solar disc (in radians).
- `denom`: Denominator factor for the CIE model, computed during initialization.
"""
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

Create a standard CIE model of sky diffuse radiance on the horizontal plane 
as described by Darula and Kittler (2002).  
See package documentation for details.

## Arguments
- `type`: The type of CIE model (1 to 15).
- `θₛ`: Zenith angle of the solar disc (in radians).
- `Φₛ`: Azimuth angle of the solar disc (in radians).
- `rtol`: Relative tolerance for the numerical integration.
- `atol`: Absolute tolerance for the numerical integration.
- `maxevals`: Maximum number of function evaluations for the numerical integration.
    
## Returns
An instance of `CIE` containing the parameters of the CIE model and the computed denominator factor.
"""
function CIE(; type = 1, θₛ = 0.0, Φₛ = 0.0, rtol = sqrt(eps(Float64)), atol = 0.0,
    maxevals = typemax(Int))
    a, b, c, d, e = CIEmodel(type)
    denom = hcubature(x -> f_denom(x[1], x[2], a, b, c, d, e, θₛ, Φₛ),
        (0.0, 0.0), (π / 2, 2π); rtol = rtol, atol = atol,
        maxevals = maxevals)
    CIE(a, b, c, d, e, θₛ, Φₛ, denom[1])
end

# Generate the parameters for the 15 CIE models
"""
    CIEmodel(type::Int)

    Get the parameters of the CIE model based on the specified type.
    There are 15 types of CIE models, each with different parameters.
    See package documentation for details.

    ## Arguments
    - `type`: An integer indicating the type of CIE model (1 to 15).

    ## Returns
    A tuple containing the parameters of the CIE model.
"""
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

# Function use to compute denom factor
"""
    f_denom(θ, Φ, a, b, c, d, e, θₛ, Φₛ)

Compute the denominator factor for the CIE model based on the zenith angle `θ`, azimuth angle `Φ`,
and the parameters of the CIE model (`a`, `b`, `c`, `d`, `e`).
See package documentation for details.

## Arguments
- `θ`: Zenith angle (in radians).
- `Φ`: Azimuth angle (in radians).
- `a`, `b`, `c`, `d`, `e`: Parameters of the CIE model.
- `θₛ`: Zenith angle of the solar disc (in radians).
- `Φₛ`: Azimuth angle of the solar disc (in radians).

## Returns
The denominator factor for the CIE model.
"""
function f_denom(θ, Φ, a, b, c, d, e, θₛ, Φₛ)
    γₛ = acos(cos(θₛ) * cos(θ) + sin(θₛ) * sin(θ) * cos(abs(Φ - Φₛ)))
    fCIE(a, b, c, d, e, θ, γₛ) * sin(θ) * cos(θ)
end

# Relative radiance with respect to any direction in the sky dome
"""
    fCIE(a, b, c, d, e, θ, γₛ)

Compute the relative radiance for a given zenith angle `θ` and solar angle `γₛ`
based on the parameters of the CIE model (`a`, `b`, `c`, `d`, `e`).
See package documentation for details.

## Arguments
- `a`, `b`, `c`, `d`, `e`: Parameters of the CIE model.
- `θ`: Zenith angle (in radians).
- `γₛ`: Solar angle (in radians).

## Returns
The relative radiance value.
"""
function fCIE(a, b, c, d, e, θ, γₛ)
    (1 + c * (exp(d * γₛ) - exp(d * π / 2)) + e * cos(γₛ)^2) * (1 + a * exp(b / cos(θ)))
end

# Radiance for a given angle normalized by horizontal irradiance
"""
    radiance(m::CIE, θ, Φ)

Calculate the radiance for a given zenith angle `θ` and azimuth angle `Φ`
based on the CIE model parameters. The radiance is normalized by the
denominator factor computed during initialization of the CIE model.
See package documentation for details.

## Arguments
- `m`: An instance of `CIE` containing the parameters of the CIE model.
- `θ`: Zenith angle (in radians).
- `Φ`: Azimuth angle (in radians).

## Returns
The radiance value normalized by the denominator factor.
"""
function radiance(m::CIE, θ, Φ)
    γₛ = acos(cos(m.θₛ) * cos(θ) + sin(m.θₛ) * sin(θ) * cos(abs(Φ - m.Φₛ)))
    fCIE(m.a, m.b, m.c, m.d, m.e, θ, γₛ) / m.denom
end

# Radiosity for a sector of the sky on horizontal plane normalized by horizontal irradiance
"""
    radiosity(m::CIE, sky::SkySectors, Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0);
    rtol = sqrt(eps(Float64)), atol = 0.0, maxevals = typemax(Int)) where {nw}

Calculate the radiosity of each section of `sky` on the horizontal plane given diffuse
irradiance on the horizontal plane (`Idif` with `nw` wavebands) assuming a CIE model.
See package documentation for details.

## Arguments
- `m`: An instance of `CIE`.
- `sky`: An instance of `SkySectors` containing the zenith and azimuth angles of each sector.
- `Idif`: A vector of diffuse irradiance values on the horizontal plane with `nw` wavebands, in W/m²/sr.
- `rtol`: Relative tolerance for the numerical integration.
- `atol`: Absolute tolerance for the numerical integration.
- `maxevals`: Maximum number of function evaluations for the numerical integration.

## Returns
An instance of `SkyDome` containing the sectors and their associated radiosity values.
"""
function radiosity(m::CIE, sky::SkySectors,
    Idif::SVector{nw, Float64} = SVector{1, Float64}(1.0);
    rtol = sqrt(eps(Float64)), atol = 0.0, maxevals = typemax(Int)) where {nw}
    I = SVector{nw, Float64}[(t = hcubature(x -> radiance(m, x[1], x[2]) * cos(x[1]) *
                                                 sin(x[1]),
        (s.θₗ, s.Φₗ), (s.θᵤ, s.Φᵤ); rtol = rtol, atol = atol,
        maxevals = maxevals)[1];
    SVector{nw, Float64}(t .* Idif[i] for i in 1:nw)) for s in sky]
    SkyDome(sky, I)
end
