### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 9eb5261b-0172-496e-8b37-af69bda54cc8
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ d5fd0104-6162-46b8-b6fb-100020af7f75
using Unitful, UnitfulAstro

# ╔═╡ 97da13b3-d291-49ee-a268-875a3cdecada
using Unitful: Hz

# ╔═╡ dd5cb545-3048-4993-bd24-9d38c88be075
using Unitful: cm

# ╔═╡ 7e05b5f6-d169-41df-bc6b-87af579cbe6e
using UnitfulAstro: foe

# ╔═╡ cb42ae83-105f-4179-b136-050494455469
using UnitfulAstro: Mpc

# ╔═╡ 97dabced-9b4f-48a6-83ca-3335b7cf0244
using CairoMakie

# ╔═╡ 351f346b-c485-444c-89a9-98ce2b2e61da
# k = 0, 2

# ╔═╡ 2d509ec6-9a55-11ef-32b6-89482b901188
#bs = [1, 2, 3, 7, 9, 10, 11]

# ╔═╡ 847d415c-a451-425f-ba7e-ee0efa2610dd
νs = logrange(1e6, 2.418e26, length = 100)Hz

# ╔═╡ a5c2a77b-c7db-4557-81d7-8ac806ebef9a
p = 2.23; # spectral index of electron distribution

# ╔═╡ 30bc52fb-f137-4233-bd4c-670cd0b58801
z = 1; # cosmological red shift

# ╔═╡ 560fb21a-d217-4239-9ae2-d8407206db51
ε_e = 0.1; # fraction of total energy density in electrons

# ╔═╡ e2dccbdd-c399-476a-9e67-6ed871c8d67e
ε_e_bar = (ε_e * (p - 2)) / (p - 1)

# ╔═╡ df809df0-0c87-45d0-8425-8e63a16afa15
ε_B = 0.01 # fraction of total energy in magnetic field

# ╔═╡ 385f3c85-2e23-4a01-bfba-d2fdf3f7b4a4
ε = 1 - (ε_e + ε_B)

# ╔═╡ 22b7badb-3c4b-425a-a4ff-4faf5896a989
n_0 = 1#/cm^3 # ambient particle density

# ╔═╡ f98dfbee-f09e-4987-83d2-60cc9f51aaeb
E = 10foe

# ╔═╡ 624292cc-9705-4654-8d3b-1ad1b6b88d7f
E_52 = E / 10foe

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
d = 6787.5Mpc |> cm

# ╔═╡ 75403e1f-db0a-4c17-ade2-69870f55f7d0
d_L28 = d / 1e28cm

# ╔═╡ 427af236-899a-4772-8a12-e8f8320e34d2
function breakcase(b, t, ν)
    # runtime dispatch instead of compile-time dispatch
    if b in (1, 2, 3, 7, 9, 10, 11)
        return breakcase(Val(b), t, ν)
    else
        error("Invalid value for b")
    end
end

# ╔═╡ 618fd47f-5c56-4914-b92a-66b14b64fe5a
function breakcase(b::Val{1}, t, ν)
    k = 0
    β_1 = 2
    β_2 = 1//3
    s = 1.64
    # ν_b = ν_sa
    # ν_sa is the self absorption frequency

    # Dummy variables for calculating ν_b
    var1 = (p-1)^(3//5) / (3p+2)^(3//5)
    var2 = 1/(1+z)
    var3 = 1/ε_e_bar
    var4 = ε_B^(1//5) * n_0^(3//5) * E_52^(1//5)

    ν_1 = 1.24e9Hz * var1 * var2 * var3 * var4

    varext1 = (p-1)^(6//5) / ((3p-1) * (3p+2)^(1//5))
    varext2 = √(1+z)
    varext3 = 1/ε_e_bar
    varext4 = ε_B^(2//5) * n_0^(7//10) * E_52^(9//10)
    varext5 = √t
    varext6 = 1/d_L28^2

    ν_1_ext = 0.647Hz * varext1 * varext2 * varext3 * varext4 * varext5 * varext6

    F_ν_1 = ν_1_ext / ((ν/ν_1).^(-s*β_1) + (ν/ν_1).^(-s*β_2))^(1/s)

    return ustrip(Hz, F_ν_1)
end

# ╔═╡ e41c8512-f0f0-4730-a337-00247adf65a6
function breakcase(b::Val{2}, t, ν)
    k = 0
    β_1 = 1//3
    β_2 = (1-p)/2
    s = 1.84 - 0.4p
    # ν_b == ν_m
    # ν_b being the peak frequency, thus ν_2 should be larger than ν_1 and ν_3

    var1 = (p - 0.67) * 1e15
    var2 = √(1+z)
    var3 = √E_52
    var4 = ε^2
    var5 = √ε_B
    var6 = t^(-3/2)

    ν_2 = 3.73Hz * var1 * var2 * var3 * var4 * var5 * var6

    F̃_2 = (1 + (ν/ν_2)^(s*(β_1 - β_2)))^(-1/s)

    return F̃_2
end

# ╔═╡ c95bd4ea-c81c-49a3-a13d-b5ea65b3eaa9
function breakcase(b::Val{3}, t, ν)
    k = 0
    β_1 = (1-p)/2
    β_2 = -p/2
    s = 1.15 - 0.06p

    var1 = (p - 0.46) * 1e13
    var2 = exp(-1.16p)
    var3 = 1/√(1+z)
    var4 = ε_B^(-3/2)
    var5 = 1/n_0
    var6 = 1/√E_52
    var7 = 1/√t

    ν_3 = 6.37Hz * var1 * var2 * var3 * var4 * var5 * var6 * var7

    F̃_3 = (1 .+ (ν/ν_3).^(s*(β_1-β_2))).^(-1/s)

    return F̃_3
end

# ╔═╡ a2c10f53-b999-42c9-abc2-fa876db201b1
function breakcase(b::Val{7}, t, ν)
    k = 0
    β_1 = 2
    β_2 = 11//8
    s = 1.9 - 0.04p

    var1 = ((3p-1)/(3p+2))^(8//5)
    var2 = 1/(1+z)^(13//10)
    var3 = ε_e_bar^(-8//5)
    var4 = ε_B^(-2//5)
    var5 = (n_0^3 * t^3 / E_52)^(1//10)

    ν_7 = 1.12e8Hz * var1 * var2 * var3 * var4 * var5

    varext1 = ((3p-1)/(3p+2))^(11//5)
    varext2 = 1/(1+z)^(1//10)
    varext3 = (ε_e_bar * ε_B)^(-4//5)
    varext4 = (n_0 * E_52^3 * t^11)^(1//10)
    varext5 = 1 / d_L28^2

    ν_7_ext = 5.27e-3Hz * varext1 * varext2 * varext3 * varext4 * varext5

    F_ν_7 = ν_7_ext ./ .√((ν/ν_7).^(-s*β_1) .* (ν/ν_7)^(-s*β_2))

    return ustrip(Hz, F_ν_7)
end

# ╔═╡ 3d66b17b-fc2a-41d8-a76a-2e40bfb8f020
function breakcase(b::Val{9}, t, ν)
    k = 0
    β_1 = -1//2
    β_2 = -p/2
    s = 3.34 - 0.82p

    var1 = p - 0.74
    var2 = √(1+z)
    var3 = ε_e_bar^2
    var4 = √ε_B
    var5 = √E_52
    var6 = t^(-3//2)

    ν_9 = 3.94e15Hz * var1 * var2 * var3 * var4 * var5 * var6

    F̃_9 = (1 .+ (ν/ν_9)^(s*(β_1 - β_2)))^(-1/s)

    return F̃_9
end

# ╔═╡ 3dc43aca-7621-4905-9df7-f43a88221016
function breakcase(b::Val{10}, t, ν)
    k = 0
    β_1 = 11//8
    β_2 = 1//3
    s = 1.213

    var1 = 1/√(1+z)
    var2 = ε_B^(6//5)
    var3 = 1//n_0
    var4 = E_52^(7//10)
    var5 = t^(-1/2)

    ν_10 = 1.32e10Hz * var1 * var2 * var3 * var4 * var5

    F̃_10 = (1 .+ (ν/ν_10).^(s*(β_1 - β_2))).^(-1/s)

    return F̃_10
end

# ╔═╡ dfa2cebc-7aaf-4539-9f13-07ce060e606e
function breakcase(b::Val{11}, t, ν)
    k = 0
    β_1 = 1//3
    β_2 = -1//2
    s = 0.597

    var1 = 1/√(1+z)
    var2 = ε_B^(-3//2)
    var3 = 1/n_0
    var4 = 1/√E_52
    var5 = 1/√t

    ν_11 = 5.86e12Hz * var1 * var2 * var3 * var4 * var5

    F̃_11 = (1 .+ (ν/ν_11)^(s^(β_1 - β_2)))^(-1/s)

    return F̃_11
end

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
    ax = Axis(f[1,1], xlabel = "log10(ν)", ylabel = "log10(F₅)")

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

# ╔═╡ Cell order:
# ╠═9eb5261b-0172-496e-8b37-af69bda54cc8
# ╠═351f346b-c485-444c-89a9-98ce2b2e61da
# ╠═2d509ec6-9a55-11ef-32b6-89482b901188
# ╠═d5fd0104-6162-46b8-b6fb-100020af7f75
# ╠═97da13b3-d291-49ee-a268-875a3cdecada
# ╠═847d415c-a451-425f-ba7e-ee0efa2610dd
# ╠═a5c2a77b-c7db-4557-81d7-8ac806ebef9a
# ╠═30bc52fb-f137-4233-bd4c-670cd0b58801
# ╠═560fb21a-d217-4239-9ae2-d8407206db51
# ╠═e2dccbdd-c399-476a-9e67-6ed871c8d67e
# ╠═df809df0-0c87-45d0-8425-8e63a16afa15
# ╠═385f3c85-2e23-4a01-bfba-d2fdf3f7b4a4
# ╠═dd5cb545-3048-4993-bd24-9d38c88be075
# ╠═22b7badb-3c4b-425a-a4ff-4faf5896a989
# ╠═7e05b5f6-d169-41df-bc6b-87af579cbe6e
# ╠═f98dfbee-f09e-4987-83d2-60cc9f51aaeb
# ╠═624292cc-9705-4654-8d3b-1ad1b6b88d7f
# ╠═43a5597d-acac-472b-8a6e-b803201ae62a
# ╠═db9b450a-c183-4650-919c-cba8bfab9b49
# ╠═cb42ae83-105f-4179-b136-050494455469
# ╠═093c21fa-d06c-4607-ac7b-f3b1e8c04139
# ╠═75403e1f-db0a-4c17-ade2-69870f55f7d0
# ╠═427af236-899a-4772-8a12-e8f8320e34d2
# ╠═618fd47f-5c56-4914-b92a-66b14b64fe5a
# ╠═e41c8512-f0f0-4730-a337-00247adf65a6
# ╠═c95bd4ea-c81c-49a3-a13d-b5ea65b3eaa9
# ╠═a2c10f53-b999-42c9-abc2-fa876db201b1
# ╠═3d66b17b-fc2a-41d8-a76a-2e40bfb8f020
# ╠═3dc43aca-7621-4905-9df7-f43a88221016
# ╠═dfa2cebc-7aaf-4539-9f13-07ce060e606e
# ╠═97dabced-9b4f-48a6-83ca-3335b7cf0244
# ╠═cb1d558d-b7ed-49fd-b64d-36a631b871e5
# ╠═a5cb1e26-8996-49c5-bf03-05724c330ec1
# ╠═c73abcf4-4e6b-4fae-bd2e-b34a2d3de4f4
