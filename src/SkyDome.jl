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
      (θₗ = s.θₗ[i], θᵤ = s.θᵤ[i],
       Φₗ = s.Φₗ[i], Φᵤ = s.Φᵤ[i])
end
length(s::SkySectors) = length(s.θₗ)
lastindex(s::SkySectors) = length(s)
iterate(s::SkySectors, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

"""
    polar(s::SkySectors)

Draw a polar graph of the sky and its discretization into sectors. Each dot 
represents the center of each sector.
"""
function polar(s::SkySectors)
    
    # Convert to data frame with ASCII names and hexadecimal degrees
    df = DataFrame(theta_lower = ustrip.(degrees.(s.θₗ)), 
                   theta_upper = ustrip.(degrees.(s.θᵤ)),
                   theta = ustrip.(degrees.(s.θₗ .+ s.θᵤ)./2),
                   phi_lower = ustrip.(degrees.(s.Φₗ)),
                   phi_upper = ustrip.(degrees.(s.Φᵤ)),
                   phi = ustrip.(degrees.(s.Φₗ .+ s.Φᵤ)./2))

    # Make sure ggplot2 is loaded in the R session
    R"if(!require(ggplot2)) install.packages('ggplot2')"
    R"library(ggplot2)"
    @rput df

    # Polar plot with ggplot2
    R"ggplot(data = df) +
        geom_rect(aes(xmin = theta_lower, xmax = theta_upper, 
                    ymin = phi_lower, ymax = phi_upper), colour = 'black', fill = 'white') +
        geom_point(aes(x = theta, y = phi)) + 
        scale_y_continuous(limits = c(0,360), expand = c(0,0),
                        breaks = c(0,45,90,135,180,225,270,315)) +
        scale_x_continuous(limits = c(0,90), expand = c(0,0),
                        breaks = c(0,15,30,45,60,75,90)) +
        coord_polar(theta = 'y') + 
        labs(x = '', y = '') +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank())"

end



################################################################################
########################### Discretization methods #############################
################################################################################

"""
    equal_angle_intervals(ntheta, nphi)

Discretize the sky into `ntheta` zenith rings of `nphi` sectors each assuming 
the same angle intervals for each sector (Δθ = π/2/ntheta and ΔΦ = 2π/nphi). 
Returns an object of type `SkySectors`. See package documentation for details.
"""
function equal_angle_intervals(ntheta, nphi)
    Δθ = π/2/ntheta*1.0rad
    ΔΦ = 2π/nphi*1.0rad
    θₗ = repeat(0.0rad:Δθ:(π/2)rad - Δθ, inner = nphi)
    θᵤ = repeat(Δθ:Δθ:(π/2)rad, inner = nphi)
    Φₗ = repeat(0.0rad:ΔΦ:(2π)rad - ΔΦ, outer = ntheta)
    Φᵤ = repeat(ΔΦ:ΔΦ:(2π)rad, outer = ntheta)

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
    Δθ = (π/2/ntheta)rad
    uθₗ = 0.0rad:Δθ:(π/2)rad - Δθ
    uθᵤ = Δθ:Δθ:(π/2)rad

    # Calculate number of azimuth sectors and ΔΦ per zenith ring 
    fac = cos.(uθₗ) - cos.(uθᵤ)
    n = ntheta*nphi
    nphis = round.(fac./sum(fac).*n)
    sum(nphis) != n && (nphis[end] = nphis[end] + n - sum(nphis))
    ΔΦs = (2π./nphis)rad

    # Generate coordinates of all sectors
    c = 1
    θₗ  = Vector{typeof(Δθ)}(undef, n)
    θᵤ = Vector{typeof(Δθ)}(undef, n)
    Φₗ  = Vector{typeof(Δθ)}(undef, n)
    Φᵤ = Vector{typeof(Δθ)}(undef, n)
    for i in 1:ntheta
        ΔΦ = ΔΦs[i]
        for j in 1:nphis[i]
            θₗ[c]  = Δθ*(i - 1)
            θᵤ[c] = Δθ*i
            Φₗ[c]  = ΔΦ*(j - 1)
            Φᵤ[c] = ΔΦ*j
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
     Φₗ = s.sectors.Φₗ[i], Φᵤ = s.sectors.Φᵤ[i],
     I = s.I[i])
end
length(s::SkyDome) = length(s.sectors.θₗ)
lastindex(s::SkyDome) = length(s)
iterate(s::SkyDome, state = 1) = state > length(s) ? nothing : (s[state], state + 1)

"""
    polar(s::SkyDome)

Draw a polar graph of the sky and its discretization into sectors with a colour
scale representing radiosity of sector as included in the `SkyDome` object.
"""
function polar(s::SkyDome)
  
  # Convert to data frame with ASCII names and hexadecimal degrees
  df = DataFrame(theta_lower = ustrip.(degrees.(s.sectors.θₗ)), 
                 theta_upper = ustrip.(degrees.(s.sectors.θᵤ)),
                 theta = ustrip.(degrees.(s.sectors.θₗ .+ s.sectors.θᵤ)./2),
                 phi_lower = ustrip.(degrees.(s.sectors.Φₗ)),
                 phi_upper = ustrip.(degrees.(s.sectors.Φᵤ)),
                 phi = ustrip.(degrees.(s.sectors.Φₗ .+ s.sectors.Φᵤ)./2),
                 I = s.I)

  # Make sure ggplot2 is loaded in the R session
  R"if(!require(ggplot2)) install.packages('ggplot2')"
  R"library(ggplot2)"
  @rput df

  # Polar plot with ggplot2
  R"ggplot(data = df) +
      geom_rect(aes(xmin = theta_lower, xmax = theta_upper, 
                  ymin = phi_lower, ymax = phi_upper, fill = I), 
                  colour = 'black') +
      geom_point(aes(x = theta, y = phi), size = 0.5) + 
      scale_y_continuous(limits = c(0,360), expand = c(0,0),
                      breaks = c(0,45,90,135,180,225,270,315)) +
      scale_x_continuous(limits = c(0,90), expand = c(0,0),
                      breaks = c(0,15,30,45,60,75,90)) +
      scale_fill_gradientn(colours = rainbow(3, rev = TRUE, alpha = 0.8)) +
      coord_polar(theta = 'y') + 
      labs(x = '', y = '') +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())"
end


"""
    radiance(m, θ, Φ)

Calculate the radiance per unit of solid angle at zenith (θ) and azimuth (Φ) 
angle normalized by diffuse radiance on the horizontal plane for a given
model of diffuse radiation `m`. See package documentation for details.
"""
function radiance(m, θ, Φ)
    nothing
end

"""
radiosity(m, sky::SkySectors, cosine = false)

Calculate the radiosity of each section of `sky`(multiplied by cos(θ) if `cosine = true`) 
normalized by diffuse radiance on the horizontal plane 
for a given model of radiation `m`. See package documentation for details.
"""
function radiosity(m, sky::SkySectors, cosine = false)
    nothing
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

# Radiance for a given angle normalized by horizontal irradiance
function radiance(m::UniformSky, θ::UN.Quantity, Φ::UN.Quantity)
    1/π/rad^2
end

# Radiosity for a sector of the sky normalized by horizontal irradiance
function radiosity(m::UniformSky, sky::SkySectors, cosine = false)
    if cosine
        I = [(s.Φᵤ - s.Φₗ)*(cos(s.θₗ)^2 - cos(s.θᵤ)^2)/2/π/rad for s in sky]
    else
        I = [(s.Φᵤ - s.Φₗ)*(cos(s.θₗ) - cos(s.θᵤ))/π/rad for s in sky]
    end
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

# Radiance for a given angle normalized by horizontal irradiance
function radiance(m::StandardSky, θ::UN.Quantity, Φ::UN.Quantity)
    (3 + 6*cos(θ))/7/π/rad^2
end

# Radiosity for a sector of the sky normalized by horizontal irradiance
function radiosity(m::StandardSky, sky::SkySectors, cosine = false)
    if cosine
        I = [(s.Φᵤ - s.Φₗ)*((4cos(s.θₗ) + 3)*cos(s.θₗ)^2 - (4cos(s.θᵤ) + 3)*cos(s.θᵤ)^2)/14/π/rad for s in sky]
    else
        I = [3*(s.Φᵤ - s.Φₗ)*((cos(s.θₗ) + 1)*cos(s.θₗ) - (cos(s.θᵤ) + 1)*cos(s.θᵤ))/7/π/rad for s in sky]
    end
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
    CIE(type::Int; θₛ = 0.0, Φₛ = 0.0, rtol = sqrt(eps(Float64)), atol = 0.0, 
maxevals = typemax(Int))

Create a standard CIE model of sky diffuse radiance as described by 
Darula and Kittler (2002). The argument `type` can have values from 1 to 15
representing the 15 standard CIE models. θₛ and Φₛ are the zenith and azimuth
angles of the solar disc. `rtol` and `atol` and `maxevals` are the relative
tolerance, absolute tolerance and maximum number of function evaluation of the
numerical integration algorithm. See package documentation for details.
"""
function CIE(type::Int; θₛ = 0.0, Φₛ = 0.0, rtol = sqrt(eps(Float64)), atol = 0.0, 
             maxevals = typemax(Int))
    a, b, c, d, e = CIEmodel(type)
    denom = hcubature(x -> f_denom(x[1]*rad, x[2]*rad, a, b, c, d, e, θₛ, Φₛ),
                      (0.0, 0.0), (π/2, 2π); rtol = rtol, atol = atol, 
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

# Function use to compute denom factor
function f_denom(θ, Φ, a, b, c, d, e, θₛ, Φₛ)
    γₛ = acos(cos(θₛ)*cos(θ) + sin(θₛ)*sin(θ)*cos(abs(Φ - Φₛ)))
    fCIE(a, b, c, d, e, θ, γₛ)*sin(θ)*cos(θ)
end

# Relative radiance with respect to any direction in the sky dome
fCIE(a, b, c, d, e, θ, γₛ) = 
             (1 + c*(exp(d*γₛ) - exp(d*π/2)) + e*cos(γₛ)^2)*(1 + a*exp(b/cos(θ)))


# Radiance for a given angle normalized by horizontal irradiance
function radiance(m::CIE, θ, Φ)
    γₛ = acos(cos(m.θₛ)*cos(θ) + sin(m.θₛ)*sin(θ)*cos(abs(Φ - m.Φₛ)))
    fCIE(m.a, m.b, m.c, m.d, m.e, θ, γₛ)/m.denom/rad^2
end

# Radiosity for a sector of the sky normalized by horizontal irradiance
function radiosity(m::CIE, sky::SkySectors, cosine = false;
                    rtol=sqrt(eps(Float64)), atol=0.0, maxevals=typemax(Int))
    if cosine
        I = [hcubature(x -> radiance(m, x[1]*rad, x[2]*rad)*cos(x[1]*rad)*sin(x[1]*rad)*rad^2,
                      (ustrip(s.θₗ),  ustrip(s.Φₗ)), 
                      (ustrip(s.θᵤ), ustrip(s.Φᵤ)); 
                      rtol = rtol, atol = atol, 
                      maxevals = maxevals)[1] for s in sky]
    else
        I = [hcubature(x -> radiance(m, x[1]*rad, x[2]*rad)*sin(x[1]*rad)*rad^2,
                      (ustrip(s.θₗ),  ustrip(s.Φₗ)), 
                      (ustrip(s.θᵤ), ustrip(s.Φᵤ)); 
                      rtol = rtol, atol = atol, 
                      maxevals = maxevals)[1] for s in sky]
    end
    SkyDome(sky, I)
end
