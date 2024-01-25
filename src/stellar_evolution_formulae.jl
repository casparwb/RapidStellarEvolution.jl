using Unitful, UnitfulAstro

# """ 
#     stellar_spin(m::T, R::T)

# Return the stellar rotation of a star with mass 'm [M⊙]' and radius 'R [R⊙]', as
# described by Hurley, Pols, & Tout 2000, eq 107-108.
# """
# function stellar_spin(m::Quantity{<:Real, mS}, R::Quantity{<:Real, RS}) where {mS, RS}
#     stellar_spin(u"Msun"(m).val, u"Rsun"(R).val)u"1/yr"
# end

# function stellar_spin(m::T, R::T) where T <: Real
#     vᵣₒₜ = 330m^3.3/(15 + m^3.45)
#     Ω = (45.35vᵣₒₜ/R)
# end


"""
hook_mass(ζ)

Initial mass above which a hook appears in the main-sequence.
Hurley et al. 2000 Equation 1
"""
function hook_mass(ζ)
    return 1.0185 + 0.16015*ζ + 0.0892*ζ^2
end

"""
hef_mass(ζ)

Return the HeF mass: maximum initial mass for which He ignites degenerately in a helium flash.
Hurley et al. 2000 Equation 2
"""
function hef_mass(ζ)
    return 1.995 + 0.25*ζ + 0.087*ζ^2
end

"""
fgb_mass(ζ)

Return the FGB mass: the maximum initial mass for which He ignites on the first giant branch.
Hurley et al. 2000 Equation 3
"""
function fbg_mass(Z)
    return (13.048*(Z/0.02)^0.06)/(1 + 0.0012*(0.02/Z)^1.27)
end

"""
bgb_lifetime(ζ)

Return t_BGB: the lifetimes to the BGB (base giant branch).
Hurley et al. 2000 Equation 4
"""
function bgb_lifetime(M, a₁, a₂, a₃, a₄, a₅)
    return (a₁ + a₂*M^4 + a₃*M^5.5 + M⁷)/(a₄*M^2 + a₅*M⁷)
end


function hook_time(M, t_BGB, a₆, a₇, a₈, a₉, a₁₀)
    μ = max(0.5, 1.0 - 0.01*max(a₆/M^a₇, a₈ + a₉/M^a₁₀))
    return μ*tBGB
end

"""
main_sequence_lifetime(ζ, t_hook, t_BGB)

Return the main sequence lifetime.
"""
function main_sequence_lifetime(ζ, t_hook, t_BGB)
    x = max(0.95, min(0.95 - 0.03*(ζ + 0.30103), 0.99))
    return max(t_hook, x*t_BGB)
end

"""
bgb_luminosity(M, a₁₁, a₁₂, a₁₃, a₁₄, a₁₅, a₁₆)

Return the luminosity at the end of the main sequence.
"""
function tams_luminosity(M, a₁₁, a₁₂, a₁₃, a₁₄, a₁₅, a₁₆=7.2)
    return (a₁₁*M^3 + a₁₂*M^4 + a₁₃*M^(a₁₆ + 1.8))/(a₁₄ + a₁₅*M^5 + M^a₁₆)
end

"""
tams_radius(M, a₁₁, a₁₂, a₁₃, a₁₄, a₁₅, a₁₆)

Return the radius at the end of the main sequence.
"""
function tams_radius(M, a₁₈, a₁₉, a₂₀, a₂₁=1.47, a₂₂=3.07, a₂₃, a₂₄, a₂₅, a₂₆=5.50)
    a₁₇ = 1.4
    M_star = a₁₇ + 0.1
    if M < a₁₇
        R_TMS = (a₁₈ + a₁₉*M^a₂₁)/(a₂₀ + M^a₂₂)
    elseif M >= M_star
        c₁ = -8.672073e-2
        R_TMS =  (c₁*M^3 + a₂₃*M^a₂₆ + a₂₄*M^(a₂₆ + 1.5))/(a₂₅ + M^5)
    end

    if M < 0.5
        return max(R_TMS, 1.5*R_TMS)
    end
    
    return R_TMS
end

"""
bgb_luminosity(M, a₁₁, a₁₂, a₁₃, a₁₄, a₁₅, a₁₆)

Return the luminosity at the base of the giant branch.
"""
function bgb_luminosity(M, a₂₇, a₂₈, a₂₉, a₃₀, a₃₁=4.60, a₃₂=6.68)
    c₂ = 9.301992
    c₃ = 4.63745
    return (a₂₇*M^a₃₁ + a₂₈*M^c₂)/(a₂₉ + a₃₀*M^c₃ + M^a₃₂)
end


function fractional_timescale(t, tMS)
    return t/tMS
end

function ms_luminosity()
    ϵ = 0.01

    η = ifelse(Z <= 0.0009, ifelse(M <= 1.1, 20, 10), 10)
end

function luminosity_perturbation(M, M_hook, a₃₃, a₃₅=0.4, a₃₆, a₃₇=0.6)

    if M <= M_hook
        return 0.0
    elseif (M_hook < M < a₃₃)
        B = luminosity_perturbation(a₃₃, M_hook, a₃₃, a₃₅, a₃₆, a₃₇)
        return B*((M - M_hook)/a₃₃ - M_hool)^0.4
    elseif M >= a₃₃
        return min(a₃₄/M^(a₃₅), a₃₆/M^a₃₇)
    end
    

end

function radius_perturbation(M, M_hook, a₃₈, a₃₉, a₄₀, a₄₁=3.57, a₄₂, a₄₃, a₄₄=1.0)

    if M <= M_hook
        return 0.0
    elseif (M_hook < M < a₄₂)
        return a₄₃*((M - M_hook)/(a₄₂ - M_hook))^0.5
    elseif (a₄₂ < M < 2.0)
        B = radius_perturbation(2, M_hook, a₃₈, a₃₉, a₄₀, a₄₁, a₄₂, a₄₃, a₄₄)
        return a₄₃ + (B - a₄₃)*((M - a₄₂)/(2 - a₄₂))^a₄₄
    elseif M >= 2
        return (a₃₈ + a₃₉*M^3.5)/(a₄₀*M^3 + M^a₄₁) - 1
    end

end

function envelope_structure(mass::Real, radius, core_mass, core_radius, stellar_type, age, Z=0.02)
    tMS, tBGB = main_sequence_lifetime(mass, Z)
    envelope_radius = convective_envelope_radius(mass, radius, core_radius, stellar_type, age, tMS, tBGB)
    envelope_mass = convective_envelope_mass(mass, core_mass, stellar_type, age, tMS, tBGB)

    return envelope_radius, envelope_mass
end

function envelope_structure(mass::Unitful.Mass, radius, core_mass, core_radius, stellar_type, age, Z=0.02)
    
    mass = ustrip(u"Msun", mass)
    radius = ustrip(u"Rsun", radius)
    core_mass = ustrip(u"Msun", core_mass)
    core_radius = ustrip(u"Rsun", core_radius)
    stellar_type = ustrip(u"stp", stellar_type)
    age = ustrip(u"Myr", age)

    R_env, M_env = envelope_structure(mass, radius, core_mass, core_radius, stellar_type, age, Z)
    return R_env*u"Rsun", M_env*u"Msun"
end


function envelope_structure(star::Particle, age, Z=0.02)
    @assert star.structure.type isa Star "Envelope structure only relevant for stars."

    envelope_structure(star.structure.m, star.structure.R, 
                       star.structure.m_core, star.structure.R_core, 
                       star.structure.type.index, age)
end


"""
Radius of a zero-age main-sequence star. From Tout et al 1996.
"""
function zero_age_main_sequence_radius(M::Real)
    θ = 1.71535900
    ι = 6.59778800
    κ = 10.08855000
    λ = 1.01249500
    μ = 0.07490166
    ν = 0.01077422
    ξ = 3.08223400
    o = 17.84778000
    Π = 0.00022582


    (θ*M^2.5 + ι*M^6.5 + κ*M^11 + λ*M^19 + μ*M^19.5)/(ν + ξ*M^2 + o*M^8.5 + M^18.5 + Π*M^19.5)
end

function zero_age_main_sequence_radius(mass::Unitful.Mass)
    zero_age_main_sequence_radius(ustrip(u"Msun", mass))
end



"""
    main_sequence_radius_035_msun(τ, Z=0.02)

Radius of a 0.35 M⊙ main-sequence star at a time τ = t/tMS.  
"""
function main_sequence_radius_035_msun(τ::Real, Z=0.02)
    M = 0.35
    ζ = log10(Z/0.02)
    ζ² = ζ^2    
    aₙ(α, β=0.0, γ=0.0, η=0.0, μ=0.0) = α + β*ζ + γ*ζ² + η*ζ^3 + μ*ζ²^2

    a₁₈′ = aₙ(2.187715e-1, -2.154437e+0, -3.768678e+0, -1.975518e+0, -3.021475e-1)
    a₁₉′ = aₙ(1.466440e+0, 1.839725e+0, 6.442199e+0, 4.023635e+0, 6.957529e-1)
    a₂₀ = aₙ(2.652091e+1, 8.178458e+1, 1.156058e+2, 7.633811e+1, 1.950698e+1)

    a₁₈ = a₁₈′*a₂₀
    a₁₉ = a₁₉′*a₂₀
    a₂₁ = aₙ(1.472103e+0, -2.947609e+0, -3.312828e+0, -9.945065e-1)
    a₂₂ = aₙ(3.071048e+0, -5.679941e+0, -9.745523e+0, -3.594543e+0)

    R_zams = zero_age_main_sequence_radius(M)
    R_tms = (a₁₈ + a₁₉*M^a₂₁)/(a₂₀ + M^a₂₂) # Hurley et al 2000 eq. 9

    Mhook = 1.0185 + 0.16015ζ + 0.0892ζ² 
    @assert Mhook >= M "Only valid for M <= Mhook right now."


    a₆₂ = aₙ(8.4300e-2, -4.7500e-2, -3.5200e-2)
    a₇₆ = aₙ(1.192334e-2, 1.083057e-2, 1.230969e+0, 1.551656e+0)
    a₇₇ = aₙ(-1.668868e-1, 5.818123e-1, -1.105027e+1, -1.668070e+1)
    a₇₈ = aₙ(7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2)
    a₇₉ = aₙ(9.409838e+0, 1.522928e+0)

    a₇₆ = max(a₇₆, -0.1015564 - 0.2161264*ζ - 0.05182516*ζ²) 
    a₇₇ = max(-0.3868776 - 0.5457078*ζ - 0.1463472*ζ², min(0.0, a₇₇))
    a₇₈ = max(0.0, min(a₇₈, 7.454 + 9.046*ζ)) 
    a₇₉ = min(a₇₉, max(2.0, -13.3 - 18.6*ζ)) 

    αR = a₆₂
    βR = 1.06

    γ = a₇₆ + a₇₇*(M - a₇₈)^a₇₉

    logRMS_over_RZAMS = αR*τ + βR * τ^10 + γ*τ^40 + 
                            (log10(R_tms/R_zams) - αR - βR - γ)*τ^3


    10^logRMS_over_RZAMS*R_zams
end


"""
    envelope_radius(mass, radius, core_radius, stellar_type)

Calculate the radius of the envelope with given mass, radius, core radius, and stellar type.
Quantities must be in units of solar mass and solar radii.
Reference Hurley et al. 2002 - DOI: 10.1046/j.1365-8711.2002.05038.x
"""
function convective_envelope_radius(mass, radius, core_radius, stellar_type, age, tMS, tBGB)

    if any(stellar_type .== (3, 5, 6, 8, 9)) # giant-like stars
        return radius - core_radius
    elseif any(stellar_type .== (1, 7))   # main sequence stars
        τ = age/tMS 
        
        R_env₀ = if mass > 1.25
                    0.0
                elseif mass < 0.35
                    radius
                else

                    R′ = main_sequence_radius_035_msun(τ)
                    # R′ is the radius of a MS star with M = 0.35 M⊙ at τ
                    return R′*sqrt(1.25 - mass)/0.9
                end

        return R_env₀*(1 - τ)^0.25
    elseif any(stellar_type .== (2, 8)) # Hertzsprung gap stars
        τ = (age - tMS)/(tBGB - tMS)
        return sqrt(τ)*(radius - core_radius)
    end

end

"""
convective_envelope_mass(mass, radius, core_radius, stellar_type)

Calculate the mass of the envelope with given stellar mass, core mass, stellar age, 
stellar main sequence lifetime, stellar base giant branch (BHG) lifetime and stellar type.
Quantities must be in units of solar mass and solar radii.
Reference Hurley et al. 2000 - https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H
"""
function convective_envelope_mass(mass, core_mass, stellar_type, age, tMS, tBGB)
    @assert stellar_types[stellar_type] isa Star "Only stars have envelopes."

    if any(stellar_type .== (1, 7)) 
        M_env₀ = if mass < 0.35
                     mass
                 elseif mass > 1.25
                     0.0
                 else
                   ( 0.35*((1.25 - mass)/0.9)^2 )
                 end
        
        τ = age/tMS
        return M_env₀*(1 - τ)^0.25
    elseif any(stellar_type .== (2, 8))
        τ = (age - tMS)/(tBGB - tMS)
        return τ*(mass - core_mass)
    else 
        return mass - core_mass
    end
end 

"""

main_sequence_lifetime(M::Real, Z)

Return the main sequence lifetime of a star with mass M [M⊙] in Myr.
Reference Hurley et al. 2000 - https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H
"""
function main_sequence_lifetime(M::Real, Z=0.02)

    ζ = log10(Z/0.02) # Hurley et al 2000 page 5

    aₙ(α, β=0.0, γ=0.0, η=0.0, μ=0.0) = α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4

    a₁ = aₙ(1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2)
    a₂ = aₙ(2.706708e3, 1.483131e3, 5.772723e2, 7.411230)
    a₃ = aₙ(1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1)
    a₄ = aₙ(4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3)
    a₅ = aₙ(3.426349e-1)
    a₆ = aₙ(1.949814e1, 1.758178, -6.008212, -4.470533)
    a₇ = aₙ(4.903830)
    a₈ = aₙ(5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3)
    a₉ = aₙ(1.312179, -3.294936e-1, 9.231860e-2, 2.610989e-2)
    a₁₀ = aₙ(8.073972e-1)

    M⁷ = M^7
    μ = max(0.5, 1.0 - 0.01*max(a₆/M^a₇, a₈ + a₉/M^a₁₀))
    x = max(0.95, min(0.95 - 0.03*(ζ + 0.30103), 0.99))

    tBGB = (a₁ + a₂*M^4 + a₃*M^5.5 + M⁷)/(a₄*M^2 + a₅*M⁷)
    t_hook = μ*tBGB

    tMS = max(t_hook, x*tBGB)

    return tMS, tBGB
end

main_sequence_lifetime(M::Unitful.Quantity, Z=0.02) = main_sequence_lifetime(ustrip(u"Msun", M), Z)
