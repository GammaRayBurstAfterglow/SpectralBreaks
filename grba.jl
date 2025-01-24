### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 9eb5261b-0172-496e-8b37-af69bda54cc8
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ bd1923a0-88e2-4817-bf58-ef4a0ebe3acf
using PlutoUI

# ╔═╡ d5fd0104-6162-46b8-b6fb-100020af7f75
begin
    using Unitful, UnitfulAstro
    using Unitful: Hz, GHz, cm
    using UnitfulAstro: foe, Mpc, mJy
end

# ╔═╡ 97dabced-9b4f-48a6-83ca-3335b7cf0244
using CairoMakie

# ╔═╡ 351f346b-c485-444c-89a9-98ce2b2e61da
# k = 0 # ISM
# k = 2 # preexplosion wind

# ╔═╡ 2d509ec6-9a55-11ef-32b6-89482b901188
#bs = [1, 2, 3, 7, 9, 10, 11] # break indices

# ╔═╡ 3a706715-c616-4ae6-ac1c-b3e85929c1ca
TableOfContents()

# ╔═╡ 847d415c-a451-425f-ba7e-ee0efa2610dd
const νs = logrange(1e6, 2.418e26, length = 100) * Hz

# ╔═╡ a5c2a77b-c7db-4557-81d7-8ac806ebef9a
const p = 2.23; # spectral index of electron distribution

# ╔═╡ 30bc52fb-f137-4233-bd4c-670cd0b58801
const z = 1; # cosmological red shift

# ╔═╡ 560fb21a-d217-4239-9ae2-d8407206db51
const εₑ = 0.1; # fraction of total energy density in electrons

# ╔═╡ e2dccbdd-c399-476a-9e67-6ed871c8d67e
const ε̄ₑ = εₑ * (p - 2) / (p - 1)

# ╔═╡ df809df0-0c87-45d0-8425-8e63a16afa15
const ε_B = 0.01 # fraction of total energy in magnetic field

# ╔═╡ 385f3c85-2e23-4a01-bfba-d2fdf3f7b4a4
const ε = 1 - (εₑ + ε_B)

# ╔═╡ 22b7badb-3c4b-425a-a4ff-4faf5896a989
const n₀ = 1#/cm^3 # ambient particle density

# ╔═╡ f98dfbee-f09e-4987-83d2-60cc9f51aaeb
const E = 10foe

# ╔═╡ 624292cc-9705-4654-8d3b-1ad1b6b88d7f
const E_52 = E / 10foe

# ╔═╡ 43a5597d-acac-472b-8a6e-b803201ae62a
#A_⋆ = []

# ╔═╡ db9b450a-c183-4650-919c-cba8bfab9b49
t_days = [
    1.157e-4,
    5.64e-2,
    5.585e-1,
    3.495,
    13.831,
] # time since gamma ray burst, in days (maybe add units?)

# ╔═╡ 093c21fa-d06c-4607-ac7b-f3b1e8c04139
d_L = 6787.5Mpc # luminosity distance

# ╔═╡ 75403e1f-db0a-4c17-ade2-69870f55f7d0
d_L28 = d_L / 1e28cm |> NoUnits

# ╔═╡ 427af236-899a-4772-8a12-e8f8320e34d2
function breakcase(b::Integer, t, ν)
    # runtime dispatch instead of compile-time dispatch
    if b in (1, 2, 3, 7, 9, 10, 11)
        return breakcase(Val(Int(b)), t, ν)
    else
        throw(ArgumentError("Invalid value for b"))
    end
end

# ╔═╡ 517fb628-3b62-4ab0-a749-72a091744ed8
md"""
Scaled frequency near the break:
```math
ϕ_b = \frac{ν}{ν_b}
```
"""

# ╔═╡ 05799234-827a-4b1c-af3f-c5eec598cbae
md"""
All assumes
```math
k = 0
```
"""

# ╔═╡ 9830bc40-d236-4629-81f4-899b308ee792
md"""
## _b_ = 1
"""

# ╔═╡ b201cfeb-bc38-4ab9-bf88-eafd9a9d5f3a
const ν₁ = let
    # Dummy variables for calculating ν_b
    var1 = (p-1)^(3//5) / (3p+2)^(3//5)
    var2 = 1/(1+z)
    var3 = 1/ε̄ₑ
    var4 = ε_B^(1//5) * n₀^(3//5) * E_52^(1//5)

    1.24e9Hz * var1 * var2 * var3 * var4
end

# ╔═╡ 618fd47f-5c56-4914-b92a-66b14b64fe5a
function breakcase(b::Val{1}, t, ν)
    k = 0
    β₁ = 2
    β₂ = 1//3
    s = 1.64
    # ν_b = ν_sa
    # ν_sa is the self absorption frequency

    varext1 = (p-1)^(6//5) / ((3p-1) * (3p+2)^(1//5))
    varext2 = √(1+z)
    varext3 = 1/ε̄ₑ
    varext4 = ε_B^(2//5) * n₀^(7//10) * E_52^(9//10)
    varext5 = √t
    varext6 = 1/d_L28^2

    ν_1_ext = 0.647mJy * varext1 * varext2 * varext3 * varext4 * varext5 * varext6

    φ₁ = ν / ν₁

    F_ν_1 = ν_1_ext / (φ₁^(-s*β₁) + φ₁^(-s*β₂))^(1/s)

    return ustrip(mJy, F_ν_1)
end

# ╔═╡ 2a69d85c-46d2-41d8-85a3-63e051b4a52e
md"""
## _b_ = 2
"""

# ╔═╡ fc817333-6a62-44a8-a6a3-5d9f90fc860f


# ╔═╡ e41c8512-f0f0-4730-a337-00247adf65a6
function breakcase(b::Val{2}, t, ν)
    k = 0
    β₁ = 1//3
    β₂ = (1-p)/2
    s = 1.84 - 0.4p
    # ν_b == ν_m
    # ν_b being the peak frequency, thus ν₂ should be larger than ν₁ and ν₃

    var1 = (p - 0.67)
    var2 = √(1+z)
    var3 = √E_52
    var4 = ε^2
    var5 = √ε_B
    var6 = t^(-3/2)

    ν₂ = 3.73e15Hz * var1 * var2 * var3 * var4 * var5 * var6

    φ₂ = ν / ν₂

    F̃_2 = (1 + φ₂^(s*(β₁ - β₂)))^(-1/s)

    return F̃_2
end

# ╔═╡ 9c6eaf0f-069d-4726-8e33-c4f34f878b80
md"""
## _b_ = 3
"""

# ╔═╡ c95bd4ea-c81c-49a3-a13d-b5ea65b3eaa9
function breakcase(b::Val{3}, t, ν)
    k = 0
    β₁ = (1-p)/2
    β₂ = -p/2
    s = 1.15 - 0.06p

    var1 = (p - 0.46) * 1e13
    var2 = exp(-1.16p)
    var3 = 1/√(1+z)
    var4 = ε_B^(-3/2)
    var5 = 1/n₀
    var6 = 1/√E_52
    var7 = 1/√t

    ν₃ = 6.37Hz * var1 * var2 * var3 * var4 * var5 * var6 * var7

    φ₃ = ν / ν₃

    F̃_3 = (1 + φ₃^(s * (β₁ - β₂)))^(-1/s)

    return F̃_3
end

# ╔═╡ 3bb51235-3d7c-4325-b463-92eab1e86e8e
md"""
## _b_ = 4
"""

# ╔═╡ 80681bd5-2101-4f4c-aa25-f3752a41e606
function breakcase(b::Val{4}, t, ν)
end

# ╔═╡ df8b59f8-507f-49fa-8825-95f6ee8aefd7
md"""
## _b_ = 7
"""

# ╔═╡ d9d18dfe-fdc3-41c1-a75c-190fb86830f4


# ╔═╡ a2c10f53-b999-42c9-abc2-fa876db201b1
function breakcase(b::Val{7}, t, ν)
    k = 0
    β₁ = 2
    β₂ = 11//8
    s = 1.9 - 0.04p

    var1 = ((3p-1)/(3p+2))^(8//5)
    var2 = 1/(1+z)^(13//10)
    var3 = ε̄ₑ^(-8//5)
    var4 = ε_B^(-2//5)
    var5 = (n₀^3 * t^3 / E_52)^(1//10)

    ν₇ = 1.12e8Hz * var1 * var2 * var3 * var4 * var5

    varext1 = ((3p-1)/(3p+2))^(11//5)
    varext2 = 1/(1+z)^(1//10)
    varext3 = (ε̄ₑ * ε_B)^(-4//5)
    varext4 = (n₀ * E_52^3 * t^11)^(1//10)
    varext5 = 1 / d_L28^2

    ν_7_ext = 5.27e-3Hz * varext1 * varext2 * varext3 * varext4 * varext5

    φ₇ = ν / ν₇

    F_ν_7 = ν_7_ext / √(φ₇^(-s*β₁) * φ₇^(-s*β₂))

    return ustrip(Hz, F_ν_7)
end

# ╔═╡ 4f7e0c99-84d4-41f2-8a0c-e76018c5cff5
md"""
## _b_ = 9
"""

# ╔═╡ 3d66b17b-fc2a-41d8-a76a-2e40bfb8f020
function breakcase(b::Val{9}, t, ν)
    k = 0
    β₁ = -1//2
    β₂ = -p/2
    s = 3.34 - 0.82p

    var1 = p - 0.74
    var2 = √(1+z)
    var3 = ε̄ₑ^2
    var4 = √ε_B
    var5 = √E_52
    var6 = t^(-3//2)

    ν_9 = 3.94e15Hz * var1 * var2 * var3 * var4 * var5 * var6

    F̃_9 = (1 + (ν/ν_9)^(s*(β₁ - β₂)))^(-1/s)

    return F̃_9
end

# ╔═╡ 260bbde7-5048-4eaa-a216-e3c10181c200
md"""
## _b_ = 10
"""

# ╔═╡ 3dc43aca-7621-4905-9df7-f43a88221016
function breakcase(b::Val{10}, t, ν)
    k = 0
    β₁ = 11//8
    β₂ = 1//3
    s = 1.213

    var1 = 1/√(1+z)
    var2 = ε_B^(6//5)
    var3 = 1//n₀
    var4 = E_52^(7//10)
    var5 = t^(-1/2)

    ν₁₀ = 1.32e10Hz * var1 * var2 * var3 * var4 * var5

    φ₁₀ = ν / ν₁₀

    F̃_10 = (1 + φ₁₀^(s * (β₁ - β₂)))^(-1/s)

    return F̃_10
end

# ╔═╡ ff4ae5ec-adc2-478d-9004-0cf18d2a6432
md"""
## _b_ = 11
"""

# ╔═╡ dfa2cebc-7aaf-4539-9f13-07ce060e606e
function breakcase(b::Val{11}, t, ν)
    k = 0
    β₁ = 1//3
    β₂ = -1//2
    s = 0.597

    var1 = 1/√(1+z)
    var2 = ε_B^(-3//2)
    var3 = 1/n₀
    var4 = 1/√E_52
    var5 = 1/√t

    ν₁₁ = 5.86e12Hz * var1 * var2 * var3 * var4 * var5

    φ₁₁ = ν / ν₁₁

    F̃_11 = (1 + φ₁₁^(s^(β₁ - β₂)))^(-1/s)

    return F̃_11
end

# ╔═╡ c02be8d3-e472-4073-9a1f-bf53f33e6558
md"""
## Plots
"""

# ╔═╡ c73abcf4-4e6b-4fae-bd2e-b34a2d3de4f4
log_F5, log_F9 = let
    F_5_big = []
    F_9_big = []
    for (i, t) in enumerate(t_days)
        @info("Computing F5, F9 for t = $t days")

        # For F5
        F_ν_1 = breakcase.(1, t, νs)
        F̃_2 = breakcase.(2, t, νs)
        F̃_3 = breakcase.(3, t, νs)

        # For F9
        F_ν_7 = breakcase.(7, t, νs)
        F̃_9 = breakcase.(9, t, νs)
        F̃_10 = breakcase.(10, t, νs)
        F̃_11 = breakcase.(11, t, νs)

        log_F5 = log10.(@. F_ν_1 * F̃_2 * F̃_3)
        log_F9 = log10.(@. F_ν_7 * F̃_9 * F̃_10 * F̃_11)
        #@debug("Calculated logged values", log_F5, log_F9)

        push!(F_5_big, log_F5)
        push!(F_9_big, log_F9)
    end
    reduce(hcat, F_5_big), reduce(hcat, F_9_big)
end

# ╔═╡ cb1d558d-b7ed-49fd-b64d-36a631b871e5
let f = Figure()
    ax = Axis(f[1,1], xlabel = "log10(ν / Hz)", ylabel = "log10(F₅)")

    for (F, t) in zip(eachcol(log_F5), t_days)
        lines!(ax, log10.(νs/Hz), F, label = "t = $t days")
    end

    axislegend(ax, position = :cb)

    f
end

# ╔═╡ a5cb1e26-8996-49c5-bf03-05724c330ec1
let f = Figure()
    ax = Axis(f[1,1], xlabel = "log10(ν)", ylabel = "log10(F₉)")

    for (F, t) in zip(eachcol(log_F9), t_days)
        lines!(ax, log10.(νs/Hz), F, label = "t = $t days")
    end

    axislegend(ax, position = :rb)

    f
end

# ╔═╡ 9dc9f11e-9a61-4509-86ad-af5cddc98f66
md"""
## Break shape function
"""

# ╔═╡ 53c6305a-860d-410b-b659-ba473d3c8403
begin
    # k = 0:
    s(::Val{0}, ::Val{1}, _::Number) = 1.64
    s(::Val{0}, ::Val{2}, p::Number) = 1.84 - 0.4p
    s(::Val{0}, ::Val{3}, p::Number) = 1.15 - 0.06p
    s(::Val{0}, ::Val{4}, p::Number) = 3.44p - 1.41
    s(::Val{0}, ::Val{5}, p::Number) = 1.47 - 0.21p
    s(::Val{0}, ::Val{6}, p::Number) = 0.94 - 0.14p
    s(::Val{0}, ::Val{7}, p::Number) = 1.99 - 0.04p
    s(::Val{0}, ::Val{8}, _::Number) = 0.907
    s(::Val{0}, ::Val{9}, p::Number) = 3.34 - 0.82p
    s(::Val{0}, ::Val{10}, _::Number) = 1.213
    s(::Val{0}, ::Val{11}, _::Number) = 0.597
    # k = 2:
    s(::Val{2}, ::Val{1}, _::Number) = 1.06
    s(::Val{2}, ::Val{2}, p::Number) = 1.76 - 0.38p
    s(::Val{2}, ::Val{3}, p::Number) = 0.8 - 0.03p
    s(::Val{2}, ::Val{4}, p::Number) = 3.63p - 1.6
    s(::Val{2}, ::Val{5}, p::Number) = 1.25 - 0.18p
    s(::Val{2}, ::Val{6}, p::Number) = 1.04 - 0.16p
    s(::Val{2}, ::Val{7}, p::Number) = 1.97 - 0.04p
    s(::Val{2}, ::Val{8}, _::Number) = 0.893
    s(::Val{2}, ::Val{9}, p::Number) = 3.68 - 0.89p
    """
        s(k, b, p)

    TODO add docstring here
    """
    s(k::Integer, b::Integer, p::Number) = s(Val(Int(k)), Val(Int(b)), p)
end

# ╔═╡ 20729c7b-2eed-4a72-b2e6-007cea58ea30
md"""
## Break frequency function
"""

# ╔═╡ 61f7cf2e-71a9-447b-9c56-8c6d09cc6cfe
breakfrequency(k, b, t) = error("To be implemented")

# ╔═╡ 33cd8279-66b5-4015-a29c-67f2aa63e120
begin
    # k = 0:
    function breakfrequency(::Val{0}, ::Val{1}, _)
        pscale = ((p-1) / (3p+2))^(3//5)
        1.24GHz * pscale * (ε_B * n₀^3 * E_52)^(1//5) / (ε̄ₑ * (1 + z))
    end

    breakfrequency(::Val{0}, ::Val{2}, t) =
        3.73e6GHz * (p - 0.67) * ε̄ₑ^2 * √(E_52 * ε_B * (1 + z) / t^3)

    breakfrequency(::Val{0}, ::Val{3}, t) =
        6.37e4GHz * (p - 0.46)*exp(-1.16p) / (n₀ * √(ε_B * (1 + z) * E_52 * t))

    breakfrequency(::Val{0}, ::Val{4}, t) =
        5.04e7GHz * (p - 1.22) * √((1+z) * ε_B * E_52 / t^3)

    function breakfrequency(::Val{0}, ::Val{5}, t)
        pscale = (4.03 - p) * exp(2.34p)
        # TODO: name these variables better
        allotherscalenumer = (ε̄ₑ^(4p-4) * ε_B^(p+2) * n₀^4 * E_52^(p+2))
        allotherscaledenom = (1+z)^(6-p) * t^(3p+2)
        3.59GHz * pscale * (allotherscalenumer/allotherscaledenom)^(1/(2p+8))
    end
    function breakfrequency(::Val{0}, ::Val{6}, t)
        scalenumer = ε̄ₑ^(4p-4) * ε_B^(p-1) * n₀^2 * E_52^(p+1)
        scaledenom = (1+z)^(7-p) * t^(3p+3)
        3.23e3GHz * (p-1.76) * (scalenumer / scaledenom)^(1/(2p+10))
    end
    function breakfrequency(::Val{0}, ::Val{7}, t)
        pscale = ((3p-1) / (3p+2))^(8//5)
        fifthscale = (ε̄ₑ^4 * ε_B)^(-2//5)
        tenthscale = ((n₀^3 * t^3) / (E_52 * (1+z)^13))^(1//10)
        0.112GHz * pscale * fifthscale * tenthscale
    end
    breakfrequency(::Val{0}, ::Val{8}, t) =
        198GHz * (n₀ * E_52)^(1//6) / √((1+z) * t)

    function breakfrequency(::Val{0}, ::Val{9}, t)
        scale = ε̄ₑ^2 * √((1+z) * ε_B * E_52 / t^3)
        3.94e6GHz * (p - 0.74) * scale
    end
    function breakfrequency(::Val{0}, ::Val{10}, t)
        scalenumer = (ε_B^12 * n₀^11 * E_52^7)^(1//10)
        scaledenom = √((1+z) * t)
        13.2GHz * scalenumer / scaledenom
    end
    breakfrequency(::Val{0}, ::Val{11}, t) =
        5.86e3GHz / (n₀ * √((1+z) * ε_B^3 * E_52 * t))
end

# ╔═╡ db423a15-7e3c-4ec8-a537-40b5ec513860
begin
    # k = 2:
    function breakfrequency(::Val{2}, ::Val{1}, t)
        pscale = ((p-1) / (3p+2))^(3//5)
        8.31GHz * pscale * (ε_B * A_star^6 / (E_52^2 * t^-3 * (1+z)^2))^(1//5)
    end
    function breakfrequency(::Val{2}, ::Val{2}, t)
        scale = √((1+z) * E_52 * ε_B / t^3) / ε̄ₑ^2
        4.02e6GHz * (p - 0.69) * scale
    end
    function breakfrequency(::Val{2}, ::Val{3}, t)
    end
    function breakfrequency(::Val{2}, ::Val{4}, t)
    end
    function breakfrequency(::Val{2}, ::Val{5}, t)
    end
    function breakfrequency(::Val{2}, ::Val{6}, t)
    end
    function breakfrequency(::Val{2}, ::Val{7}, t)
    end
    function breakfrequency(::Val{2}, ::Val{8}, t)
    end
    function breakfrequency(::Val{2}, ::Val{9}, t)
    end
end


# ╔═╡ 58e961f9-218e-4780-acec-81c7e536c20c
md"""
## Slope function
"""

# ╔═╡ 703bd92b-a962-42f3-8e78-b32bb8c60feb
begin
    slopebelow(::Val{1}, ::Number) = 2
    slopebelow(::Val{2}, ::Number) = 1//3
    slopebelow(::Val{3}, p::Number) = (1 - p) // 2
    slopebelow(::Val{4}, ::Number) = 2
    slopebelow(::Val{5}, ::Number) = 5//2
    slopebelow(::Val{6}, ::Number) = 5//2
    slopebelow(::Val{7}, ::Number) = 2
    slopebelow(::Val{8}, ::Number) = 11//8
    slopebelow(::Val{9}, ::Number) = -1//2
    slopebelow(::Val{10}, ::Number) = 11//8
    slopebelow(::Val{11}, ::Number) = 1//3
    """
        slopebelow(b, p)

    TODO write docstring here about β₁
    """
    slopebelow(b::Integer, p::Number) = slopebelow(Val(Int(b)), p)
end

# ╔═╡ f73209ba-f495-4e6a-9e07-1dcb30579f8f
begin
    slopeabove(::Val{1}, ::Number) = 1//3
    slopeabove(::Val{2}, p::Number) = (1-p)//2
    slopeabove(::Val{3}, p::Number) = -p//2
    slopeabove(::Val{4}, ::Number) = 5//2
    slopeabove(::Val{5}, p::Number) = (1-p)//2
    slopeabove(::Val{6}, p::Number) = -p//2
    slopeabove(::Val{7}, ::Number) = 11//8
    slopeabove(::Val{8}, ::Number) = -1//2
    slopeabove(::Val{9}, p::Number) = -p//2
    slopeabove(::Val{10}, ::Number) = 1//3
    slopeabove(::Val{11}, ::Number) = -1//2
    """
        slopeabove(b, p)

    TODO write docstring here about β₂
    """
    slopeabove(b::Integer, p::Number) = slopeabove(Val(Int(b)), p)
end

# ╔═╡ 9f49f5c5-8499-45ff-8eb8-abdf7838574c
md"""
## Equation reference
"""

# ╔═╡ 407bdeaa-da92-4263-8e04-ac5f113e4d15
md"""
TODO: rewrite equations from Granot and Sari here
"""

# ╔═╡ Cell order:
# ╠═9eb5261b-0172-496e-8b37-af69bda54cc8
# ╠═351f346b-c485-444c-89a9-98ce2b2e61da
# ╠═2d509ec6-9a55-11ef-32b6-89482b901188
# ╠═bd1923a0-88e2-4817-bf58-ef4a0ebe3acf
# ╠═3a706715-c616-4ae6-ac1c-b3e85929c1ca
# ╠═d5fd0104-6162-46b8-b6fb-100020af7f75
# ╠═847d415c-a451-425f-ba7e-ee0efa2610dd
# ╠═a5c2a77b-c7db-4557-81d7-8ac806ebef9a
# ╠═30bc52fb-f137-4233-bd4c-670cd0b58801
# ╠═560fb21a-d217-4239-9ae2-d8407206db51
# ╠═e2dccbdd-c399-476a-9e67-6ed871c8d67e
# ╠═df809df0-0c87-45d0-8425-8e63a16afa15
# ╠═385f3c85-2e23-4a01-bfba-d2fdf3f7b4a4
# ╠═22b7badb-3c4b-425a-a4ff-4faf5896a989
# ╠═f98dfbee-f09e-4987-83d2-60cc9f51aaeb
# ╠═624292cc-9705-4654-8d3b-1ad1b6b88d7f
# ╠═43a5597d-acac-472b-8a6e-b803201ae62a
# ╠═db9b450a-c183-4650-919c-cba8bfab9b49
# ╠═093c21fa-d06c-4607-ac7b-f3b1e8c04139
# ╠═75403e1f-db0a-4c17-ade2-69870f55f7d0
# ╠═427af236-899a-4772-8a12-e8f8320e34d2
# ╟─517fb628-3b62-4ab0-a749-72a091744ed8
# ╟─05799234-827a-4b1c-af3f-c5eec598cbae
# ╟─9830bc40-d236-4629-81f4-899b308ee792
# ╠═b201cfeb-bc38-4ab9-bf88-eafd9a9d5f3a
# ╠═618fd47f-5c56-4914-b92a-66b14b64fe5a
# ╟─2a69d85c-46d2-41d8-85a3-63e051b4a52e
# ╠═fc817333-6a62-44a8-a6a3-5d9f90fc860f
# ╠═e41c8512-f0f0-4730-a337-00247adf65a6
# ╟─9c6eaf0f-069d-4726-8e33-c4f34f878b80
# ╠═c95bd4ea-c81c-49a3-a13d-b5ea65b3eaa9
# ╟─3bb51235-3d7c-4325-b463-92eab1e86e8e
# ╠═80681bd5-2101-4f4c-aa25-f3752a41e606
# ╟─df8b59f8-507f-49fa-8825-95f6ee8aefd7
# ╠═d9d18dfe-fdc3-41c1-a75c-190fb86830f4
# ╠═a2c10f53-b999-42c9-abc2-fa876db201b1
# ╟─4f7e0c99-84d4-41f2-8a0c-e76018c5cff5
# ╠═3d66b17b-fc2a-41d8-a76a-2e40bfb8f020
# ╟─260bbde7-5048-4eaa-a216-e3c10181c200
# ╠═3dc43aca-7621-4905-9df7-f43a88221016
# ╟─ff4ae5ec-adc2-478d-9004-0cf18d2a6432
# ╠═dfa2cebc-7aaf-4539-9f13-07ce060e606e
# ╟─c02be8d3-e472-4073-9a1f-bf53f33e6558
# ╠═97dabced-9b4f-48a6-83ca-3335b7cf0244
# ╠═cb1d558d-b7ed-49fd-b64d-36a631b871e5
# ╠═a5cb1e26-8996-49c5-bf03-05724c330ec1
# ╠═c73abcf4-4e6b-4fae-bd2e-b34a2d3de4f4
# ╟─9dc9f11e-9a61-4509-86ad-af5cddc98f66
# ╠═53c6305a-860d-410b-b659-ba473d3c8403
# ╟─20729c7b-2eed-4a72-b2e6-007cea58ea30
# ╠═61f7cf2e-71a9-447b-9c56-8c6d09cc6cfe
# ╠═33cd8279-66b5-4015-a29c-67f2aa63e120
# ╠═db423a15-7e3c-4ec8-a537-40b5ec513860
# ╟─58e961f9-218e-4780-acec-81c7e536c20c
# ╠═703bd92b-a962-42f3-8e78-b32bb8c60feb
# ╠═f73209ba-f495-4e6a-9e07-1dcb30579f8f
# ╟─9f49f5c5-8499-45ff-8eb8-abdf7838574c
# ╟─407bdeaa-da92-4263-8e04-ac5f113e4d15
