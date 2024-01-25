
struct SECoefficient{T<:Real}
    α::T
    β::T
    γ::T
    η::T
    μ::T
    SECoefficient(α=0.0, β=0.0, γ=0.0, η=0.0, μ=0.0) = new{typeof(α)}(α, β, γ, η, μ)
end

function xₙ(ζ, α, β, γ, η, μ)
    return α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4
end

function xₙ(ζ, coeff::SECoefficient)
    return xₙ(ζ, coeff.α, coeff.β, coeff.γ, coeff.η, coeff.μ)
end

function get_a_coefficient(n, Z, ζ, σ, ρ)
    if n == 11
        a₁₁_prime = xₙ(ζ, a_parameters[11])
        a₁₄ = xₙ(ζ, a_parameters[14])
        return a₁₁_prime*a₁₄
    elseif n == 12
        a₁₂_prime = xₙ(ζ, a_parameters[12])
        a₁₄ = xₙ(ζ, a_parameters[14])
        return a₁₂_prime*a₁₄
    elseif n == 17
        log_a₁₇ = max(0.097 - 0.1072*(σ + 3), max(0.097, min(0.1461, 0.1461 + 0.1237*(σ + 2))))
        return 10^log_a₁₇
    elseif n == 18
        a₁₈_prime = xₙ(ζ, a_parameters[18])
        a₂₀ = xₙ(ζ, a_parameters[20])
        return a₁₈_prime*a₂₀
    elseif n == 19
        a₁₉_prime = xₙ(ζ, a_parameters[19])
        a₂₀ = xₙ(ζ, a_parameters[20])
        return a₁₉_prime*a₂₀
    elseif n == 29
        a₂₉_prime = xₙ(ζ, a_parameters[29])
        a₃₂ = xₙ(ζ, a_parameters[32])
        return a₂₉_prime^a₃₂
    elseif n == 33
        a₃₃ =  min(1.4, 1.5135 + 0.3769ζ)
        return max(0.6355 - 0.4192ζ, max(1.25, a₃₃))
    elseif n == 42
        a₄₂ = xₙ(ζ, a_parameters[42])
        return min(1.25, max(1.2, a₄₂))
    elseif n == 44
        a₄₄ = xₙ(ζ, a_parameters[44])
        return min(1.3, max(0.45, a₄₄))
    elseif n == 49
        a₄₉ = xₙ(ζ, a_parameters[49])
        return max(a₄₉, 0.145)
    elseif n == 50
        a₅₀ = xₙ(ζ, a_parameters[50])
        return min(a₅₀, 0.306 + 0.053ζ)
    elseif n == 51
        a₅₁ = xₙ(ζ, a_parameters[51])
        return min(a₅₁, 0.3625 + 0.062ζ)
    elseif n == 52
        a₅₂ = xₙ(ζ, a_parameters[52])
        a₅₂ = max(a₅₂, 0.9)
        return ifelse(Z > 0.01, min(a₅₂, 1.0), a₅₂)
    elseif n == 53
        a₅₃ = xₙ(ζ, a_parameters[53])
        a₅₃ = min(a₅₃, 1.0)
        return ifelse(Z > 0.01, min(a₅₃, 1.1), a₅₃)
    elseif n == 57
        a₅₇ = xₙ(ζ, a_parameters[57])
        a₅₇ = min(1.4, a₅₇)
        return max(0.6355 - 0.4192ζ, max(1.25, a₅₇))
    elseif n == 62
        a₆₂ = xₙ(ζ, a_parameters[62])
        return max(0.065, a₆₂)
    elseif n == 63
        a₆₃ = xₙ(ζ, a_parameters[63])
        return ifelse(Z < 0.004, min(0.055, a₆₃), a₆₃)
    elseif n == 64
        a₆₆ = get_a_coefficient(66, Z, ζ, σ, ρ)
        a₆₈ = get_a_coefficient(68, Z, ζ, σ, ρ)

        if a₆₆ > a₆₈
            a₅₈ = get_a_coefficient(58, Z, ζ, σ, ρ)
            a₅₉ = get_a_coefficient(59, Z, ζ, σ, ρ)
            a₆₀ = get_a_coefficient(60, Z, ζ, σ, ρ)
            a₆₁ = get_a_coefficient(61, Z, ζ, σ, ρ)

            M = a₆₆
            αR = (a₅₈*M^a₆₀)/(a₅₉*M^a₆₁)
            return αR
        else
            a₆₄ = xₙ(ζ, a_parameters[64])
            return max(0.091, min(0.121, a₆₄))
        end
    elseif n == 66
        a₆₆ = xₙ(ζ, a_parameters[66])
        a₆₆ = max(a₆₆, min(1.6, -0.308 - 1.046ζ))
        return max(0.8, min(0.8 - 2ζ, a₆₆))
    elseif n == 68
        a₆₆ = get_a_coefficient(66, Z, ζ, σ, ρ)

        a₆₈ = xₙ(ζ, a_parameters[68])
        a₆₈ = max(0.9, min(a₆₈, 1.0))
        return min(a₆₈, a₆₆)
    elseif n == 72
        a₇₂ = xₙ(ζ, a_parameters[72])
        a₇₂ = ifelse(Z > 0.01, max(a₇₂, 0.95), a₇₂)
    elseif n == 74
        a₇₄ = xₙ(ζ, a_parameters[74])
        return max(1.4, min(a₇₄, 1.6))
    elseif n == 75
        a₇₅ = xₙ(ζ, a_parameters[74])
        a₇₅ = max(1.0, min(a₇₅, 1.27))
        return max(a₇₅, 0.6355 - 0.4192ζ)
    elseif n == 76
        a₇₆ = xₙ(ζ, a_parameters[76])
        return max(a₇₆, -0.1015564 - 0.2161264ζ - 0.05182516ζ^2)
    elseif n == 77
        a₇₇ = xₙ(ζ, a_parameters[77])
        return max(-0.3868776 - 0.5457078ζ - 0.1463472ζ^2 , min(0.0, a₇₇))
    elseif n == 78
        a₇₈ = xₙ(ζ, a_parameters[78])
        return max(0.0, min(a₇₈, 7.454 + 9.046ζ))
    elseif n == 79
        a₇₉ = xₙ(ζ, a_parameters[79])
        return min(a₇₉, max(2.0, -13.3 - 18.6ζ))
    elseif n == 80
        a₈₀ = xₙ(ζ, a_parameters[80])
        return max(0.0585542, a₈₀)
    elseif n == 81
        a₈₁ = xₙ(ζ, a_parameters[81])
        return min(1.5, max(0.4, a₈₁))
    else
        return xₙ(ζ, a_parameters[n])
    end
end 



function get_b_coefficient(n, Z, ζ, σ, ρ)

    if n == 1
        b1 = xₙ(ζ, b_parameters[1])
        return min(0.54, b1)
    elseif n == 2
        b2 = 10^(-4.6739 - 0.9394σ)
        return min(max(b2, -0.04167 + 55.67*Z), 0.4771 - 9329.21*Z^2.94)
    elseif n == 3
        b3_prime = max(-0.1451, -2.2794, -1.5175*σ, -0.254*σ^2)
        b3 = 10^b3_prime
        return ifelse(Z > 0.004, max(b3, 0.7307 + 14265.1*Z^3.395), b3)
    elseif n == 4
        b4 = xₙ(ζ, b_parameters[4])
        return b4 + 0.1231572*ζ^5
    elseif n == 6
        b6 = xₙ(ζ, b_parameters[6])
        return b6 + 0.01640687*ζ^5
    elseif n == 11
        b11_prime = xₙ(ζ, b_parameters[11])
        return b11_prime^2
    elseif n == 13
        b13_prime = xₙ(ζ, b_parameters[13])
        return b13_prime^2
    elseif n == 14
        b14_prime = xₙ(ζ, b_parameters[14])
        b15 = get_b_coefficient(15, Z, ζ, σ, ρ)
        return b14_prime^b15
    elseif n == 16
        b16_prime = xₙ(ζ, b_parameters[16])
        b15 = get_b_coefficient(15, Z, ζ, σ, ρ)
        return b16_prime^b15
    elseif n == 17
        b17 = 1.0
        return ifelse(ζ > -1, 1 - 0.3880523*(ζ + 1)^2.862149, b17)
    elseif n == 24
        b24_prime = xₙ(ζ, b_parameters[24])
        b28 = xₙ(ζ, b_parameters[28])
        return b24_prime^b28
    elseif n == 26
        return 5 - 0.09138012*Z^-0.3671407
    elseif n == 27
        b27_prime = xₙ(ζ, b_parameters[27])
        b28 = xₙ(ζ, b_parameters[28])
        return b27_prime^(2*b28)
    elseif n == 31
        b31_prime = xₙ(ζ, b_parameters[31])
        b33 = get_b_coefficient(33, Z, ζ, σ, ρ)
        return b31_prime^b33
    elseif n == 34
        b34_prime = xₙ(ζ, b_parameters[34])
        b33 = get_b_coefficient(33, Z, ζ, σ, ρ)
        return b34_prime^b33
    elseif n == 36
        b36_prime = xₙ(ζ, b_parameters[36])
        return b36_prime^4
    elseif n == 37
        b37_prime = xₙ(ζ, b_parameters[37])
        return 4*b37_prime
    elseif n == 38
        b38_prime = xₙ(ζ, b_parameters[38])
        return b38_prime^4
    elseif n == 40
        b40 = xₙ(ζ, b_parameters[40])
        return max(b40, 1.0)
    elseif n == 41
        b41_prime = xₙ(ζ, b_parameters[41])
        b42 = get_b_coefficient(42, Z, ζ, σ, ρ)
        return b41_prime^b42
    elseif n == 44
        b44_prime = xₙ(ζ, b_parameters[44])
        return b44_prime^5
    elseif n == 45
        return ifelse(ρ <= 0.0, 1.0, 1 - (2.47162ρ - 5.401682ρ^2 + 3.247361ρ^3))
    elseif n == 46
        b46 = xₙ(ζ, b_parameters[46])
        M_HeF = 1.995 + 0.25*ζ + 0.087*ζ^2
        M_FGB = (13.048*(Z/0.02)^0.06)/(1 + 0.0012*(0.02/Z)^1.27)
        return -1.0*b46*log10(M_HeF/M_FGB)
    elseif n == 47
        return 1.127733ρ + 0.2344416ρ^2 - 0.3793726ρ^3
    elseif n == 51
        b51_prime = xₙ(ζ, b_parameters[51])
        return b51_prime - 0.1343789*ζ^5
    elseif n == 53
        b53_prime = xₙ(ζ, b_parameters[53])
        return b53_prime + 0.4426929*ζ^5
    elseif n == 55
        b55 = xₙ(ζ, b_parameters[55])
        return min(0.99164 - 743.123*Z^2.83, b55)
    elseif n == 56
        b56_prime = xₙ(ζ, b_parameters[56])
        return b56_prime + 0.1140142*ζ^5
    elseif n == 57
        b57_prime = xₙ(ζ, b_parameters[57])
        return b57_prime - 0.01308728*ζ^5
    else
        return xₙ(ζ, b_parameters[n])
    end
end

const a_parameters = Dict{Int, SECoefficient}(1 => SECoefficient(1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2),
                                                2 => SECoefficient(2.706708e3, 1.483131e3, 5.772723e2, 7.411230),
                                                3 => SECoefficient(1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1),
                                                4 => SECoefficient(4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3),
                                                5 => SECoefficient(3.426349e-1),
                                                6 => SECoefficient(1.949814e1, 1.758178, -6.008212, -4.470533),
                                                7 => SECoefficient(4.903830),
                                                8 => SECoefficient(5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3),
                                                9 => SECoefficient(1.312179, -3.294936e-1, 9.231860e-2, 2.610989e-2),
                                                10 => SECoefficient(8.073972e-1),
                                                11 => SECoefficient(1.031538e+0, -2.434480e-1, 7.732821e+0, 6.460705e+0, 1.374484e+0),
                                                12 => SECoefficient(1.043715e+0, -1.577474e+0, -5.168234e+0, -5.596506e+0, -1.299394e+0),
                                                13 => SECoefficient(7.859573e+2, -8.542048e+0, -2.642511e+1, -9.585707e+0),
                                                14 => SECoefficient(3.858911e+3, 2.459681e+3, -7.630093e+1, -3.486057e+2, -4.861703e+1),
                                                15 => SECoefficient(2.888720e+2, 2.952979e+2, 1.850341e+2, 3.797254e+1),
                                                16 => SECoefficient(7.196580e+0, 5.613746e-1, 3.805871e-1, 8.398728e-2),

                                                18 => SECoefficient(2.187715e-1, -2.154437e+0, -3.768678e+0, -1.975518e+0, -3.021475e-1),
                                                19 => SECoefficient(1.466440e+0, 1.839725e+0, 6.442199e+0, 4.023635e+0, 6.957529e-1),
                                                20 => SECoefficient(2.652091e+1, 8.178458e+1, 1.156058e+2, 7.633811e+1, 1.950698e+1),
                                                21 => SECoefficient(1.472103e+0, -2.947609e+0, -3.312828e+0, -9.945065e-1),
                                                22 => SECoefficient(3.071048e+0, -5.679941e+0, -9.745523e+0, -3.594543e+0),
                                                23 => SECoefficient(2.617890e+0, 1.019135e+0, -3.292551e-2, -7.445123e-2),
                                                24 => SECoefficient(1.075567e-2, 1.773287e-2, 9.610479e-3, 1.732469e-3),
                                                25 => SECoefficient(1.476246e+0, 1.899331e+0, 1.195010e+0, 3.035051e-1),
                                                26 => SECoefficient(5.502535e+0, -6.601663e-2, 9.968707e-2, 3.599801e-2),
                                                27 => SECoefficient(9.511033e+1, 6.819618e+1, -1.045625e+1, -1.474939e+1),
                                                28 => SECoefficient(3.113458e+1, 1.012033e+1, -4.650511e+0, -2.463185e+0),
                                                29 => SECoefficient(1.413057e+0, 4.578814e-1, -6.850581e-2, -5.588658e-2),
                                                30 => SECoefficient(3.910862e+1, 5.196646e+1, 2.264970e+1, 2.873680e+0),
                                                31 => SECoefficient(4.597479e+0, -2.855179e-1, 2.709724e-1),
                                                32 => SECoefficient(6.682518e+0, 2.827718e-1, -7.294429e-2),
                                                
                                                34 => SECoefficient(1.910302e-1, 1.158624e-1, 3.348990e-2, 2.599706e-3),
                                                35 => SECoefficient(3.931056e-1, 7.277637e-2, -1.366593e-1, -4.508946e-2),
                                                36 => SECoefficient(3.267776e-1, 1.204424e-1, 9.988332e-2, 2.455361e-2),
                                                37 => SECoefficient(5.990212e-1, 5.570264e-2, 6.207626e-2, 1.777283e-2),
                                                38 => SECoefficient(7.330122e-1, 5.192827e-1, 2.316416e-1, 8.346941e-3),
                                                39 => SECoefficient(1.172768e+0, -1.209262e-1, -1.193023e-1, -2.859837e-2),
                                                40 => SECoefficient(3.982622e-1, -2.296279e-1, -2.262539e-1, -5.219837e-2),
                                                41 => SECoefficient(3.571038e+0, -2.223625e-2, -2.611794e-2, -6.359648e-3),
                                                42 => SECoefficient(1.9848e+0, 1.1386e+0, 3.5640e-1),
                                                43 => SECoefficient(6.300e-2, 4.810e-2, 9.840e-3),
                                                44 => SECoefficient(1.200e+0, 2.450e+0),
                                                45 => SECoefficient(2.321400e-1, 1.828075e-3, -2.232007e-2, -3.378734e-3),
                                                46 => SECoefficient(1.163659e-2, 3.427682e-3, 1.421393e-3, -3.710666e-3),
                                                47 => SECoefficient(1.048020e-2, -1.231921e-2, -1.686860e-2, -4.234354e-3),
                                                48 => SECoefficient(1.555590e+0, -3.223927e-1, -5.197429e-1, -1.066441e-1),
                                                49 => SECoefficient(9.7700e-2, -2.3100e-1, -7.5300e-2),
                                                50 => SECoefficient(2.4000e-1, 1.8000e-1, 5.9500e-1),
                                                51 => SECoefficient(3.3000e-1, 1.3200e-1, 2.1800e-1),
                                                52 => SECoefficient(1.1064e+0, 4.1500e-1, 1.8000e-1),
                                                53 => SECoefficient(1.1900e+0, 3.7700e-1, 1.7600e-1),
                                                54 => SECoefficient(3.855707e-1, -6.104166e-1, 5.676742e+0, 1.060894e+1, 5.284014e+0),
                                                55 => SECoefficient(3.579064e-1, -6.442936e-1, 5.494644e+0, 1.054952e+1, 5.280991e+0),
                                                56 => SECoefficient(9.587587e-1, 8.777464e-1, 2.017321e-1),
                                                57 => SECoefficient(1.5135e+0, 3.7690e-1),
                                                58 => SECoefficient(4.907546e-1, -1.683928e-1, -3.108742e-1, -7.202918e-2),
                                                59 => SECoefficient(4.537070e+0, -4.465455e+0, -1.612690e+0, -1.623246e+0),
                                                60 => SECoefficient(1.796220e+0, 2.814020e-1, 1.423325e+0, 3.421036e-1),
                                                61 => SECoefficient(2.256216e+0, 3.773400e-1, 1.537867e+0, 4.396373e-1),
                                                62 => SECoefficient(8.4300e-2, -4.7500e-2, -3.5200e-2),
                                                63 => SECoefficient(7.3600e-2, 7.4900e-2, 4.4260e-2),
                                                64 => SECoefficient(1.3600e-1, 3.5200e-2),
                                                65 => SECoefficient(1.564231e-3, 1.653042e-3, -4.439786e-3, -4.951011e-3, -1.216530e-03),
                                                66 => SECoefficient(1.4770e+0, 2.9600e-1),
                                                67 => SECoefficient(5.210157e+0, -4.143695e+0, -2.120870e+0),
                                                68 => SECoefficient(1.1160e+0, 1.6600e-1),
                                                69 => SECoefficient(1.071489e+0, -1.164852e-1, -8.623831e-2, -1.582349e-2),
                                                70 => SECoefficient(7.108492e-1, 7.935927e-1, 3.926983e-1, 3.622146e-2),
                                                71 => SECoefficient(3.478514e+0, -2.585474e-2, -1.512955e-2, -2.833691e-3),
                                                72 => SECoefficient(9.132108e-1, -1.653695e-1, 0.0, 3.636784e-2),
                                                73 => SECoefficient(3.969331e-3, 4.539076e-3, 1.720906e-3, 1.897857e-4),
                                                74 => SECoefficient(1.600e+0, 7.640e-1, 3.322e-1),
                                                75 => SECoefficient(8.109e-1, -6.282e-1),
                                                76 => SECoefficient(1.192334e-2, 1.083057e-2, 1.230969e+0, 1.551656e+0),
                                                77 => SECoefficient(-1.668868e-1, 5.818123e-1, -1.105027e+1, -1.668070e+1),
                                                78 => SECoefficient(7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2),
                                                79 => SECoefficient(9.409838e+0, 1.522928e+0),
                                                80 => SECoefficient(-2.7110e-1, -5.7560e-1, -8.3800e-2),
                                                81 => SECoefficient(2.4930e+0, 1.1475e+0))



const b_parameters = Dict{Int, SECoefficient}(1 => SECoefficient(3.9700e-1, 2.8826e-1, 5.2930e-1),

                                                4 => SECoefficient(9.960283e-1, 8.164393e-1, 2.383830e+0, 2.223436e+0, 8.638115e-1),
                                                5 => SECoefficient(2.561062e-1, 7.072646e-2, -5.444596e-2, -5.798167e-2, -1.349129e-2),
                                                6 => SECoefficient(1.157338e+0, 1.467883e+0, 4.299661e+0, 3.130500e+0, 6.992080e-1),
                                                7 => SECoefficient(4.022765e-1, 3.050010e-1, 9.962137e-1, 7.914079e-1, 1.728098e-1),
                                                9 => SECoefficient(2.751631e+3, 3.557098e+2),
                                                10 => SECoefficient(-3.820831e-2, 5.872664e-2),
                                                11 => SECoefficient(1.071738e+2, -8.970339e+1, -3.949739e+1),
                                                12 => SECoefficient(7.348793e+2, -1.531020e+2, -3.793700e+1),
                                                13 => SECoefficient(9.219293e+0, -2.005865e+0, -5.561309e-1),
                                                14 => SECoefficient(2.917412e+0, 1.575290e+0, 5.751814e-1),
                                                15 => SECoefficient(3.629118e+0, -9.112722e-1, 1.042291e+0),
                                                16 => SECoefficient(4.916389e+0, 2.862149e+0, 7.844850e-1),

                                                18 => SECoefficient(5.496045e+1, -1.289968e+1, 6.385758e+0),
                                                19 => SECoefficient(1.832694e+0, -5.766608e-2, 5.696128e-2),
                                                20 => SECoefficient(1.211104e+2),
                                                21 => SECoefficient(2.214088e+2, 2.187113e+2, 1.170177e+1, -2.635340e+1),
                                                22 => SECoefficient(2.063983e+0, 7.363827e-1, 2.654323e-1, -6.140719e-2),
                                                23 => SECoefficient(2.003160e+0, 9.388871e-1, 9.656450e-1, 2.362266e-1),
                                                24 => SECoefficient(1.609901e+1, 7.391573e+0, 2.277010e+1, 8.334227e+0),
                                                25 => SECoefficient(1.747500e-1, 6.271202e-2, -2.324229e-2, -1.844559e-2),

                                                27 => SECoefficient(2.752869e+0, 2.729201e-2, 4.996927e-1, 2.496551e-1),
                                                28 => SECoefficient(3.518506e+0, 1.112440e+0, -4.556216e-1, -2.179426e-1),
                                                29 => SECoefficient(1.626062e+2, -1.168838e+1, -5.498343e+0),
                                                30 => SECoefficient(3.336833e-1, -1.458043e-1, -2.011751e-2),
                                                31 => SECoefficient(7.425137e+1, 1.790236e+1, 3.033910e+1, 1.018259e+1),
                                                32 => SECoefficient(9.268325e+2, -9.739859e+1, -7.702152e+1, -3.158268e+1),
                                                33 => SECoefficient(2.474401e+0, 3.892972e-1),
                                                34 => SECoefficient(1.127018e+1, 1.622158e+0, -1.443664e+0, -9.474699e-1),
                                                
                                                36 => SECoefficient(1.445216e-1, -6.180219e-2, 3.093878e-2, 1.567090e-2),
                                                37 => SECoefficient(1.304129e+0, 1.395919e-1, 4.142455e-3, -9.732503e-3),
                                                38 => SECoefficient(5.114149e-1, -1.160850e-2),
                                                39 => SECoefficient(1.314955e+2, 2.009258e+1, -5.143082e-1, -1.379140e+0),
                                                40 => SECoefficient(1.823973e+1, -3.074559e+0, -4.307878e+0),
                                                41 => SECoefficient(2.327037e+0, 2.403445e+0, 1.208407e+0, 2.087263e-1),
                                                42 => SECoefficient(1.997378e+0, -8.126205e-1),
                                                43 => SECoefficient(1.079113e-1, 1.762409e-2, 1.096601e-2, 3.058818e-3),
                                                44 => SECoefficient(2.327409e+0, 6.901582e-1, -2.158431e-1, -1.084117e-1),
                                                46 => SECoefficient(2.214315e+0, -1.975747e+0),
                                                47 => SECoefficient(5.072525e+0, 1.146189e+1, 6.961724e+0, 1.316965e+0),
                                                48 => SECoefficient(5.072525e+0, 1.146189e+1, 6.961724e+0, 1.316965e+0),
                                                49 => SECoefficient(5.139740e+0),

                                                51 => SECoefficient(1.125124e+0, 1.306486e+0, 3.622359e+0, 2.601976e+0, 3.031270e-1),
                                                52 => SECoefficient(3.349489e-1, 4.531269e-3, 1.131793e-1, 2.300156e-1, 7.632745e-2),
                                                53 => SECoefficient(1.467794e+0, 2.798142e+0, 9.455580e+0, 8.963904e+0, 3.339719e+0),
                                                54 => SECoefficient(4.658512e-1, 2.597451e-1, 9.048179e-1, 7.394505e-1, 1.607092e-1),
                                                55 => SECoefficient(1.0422e+0, 1.3156e-1, 4.5000e-2),
                                                56 => SECoefficient(1.110866e+0, 9.623856e-1, 2.735487e+0, 2.445602e+0, 8.826352e-1),
                                                57 => SECoefficient(-1.584333e-1, -1.728865e-1, -4.461431e-1, -3.925259e-1, -1.276203e-1)
                                               )
