### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ aee7b830-cd3f-11ef-2c81-6b39b7339e04
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ d96fe018-58aa-4e2f-b6d5-272f66c2ce8b
using StaticArrays

# ╔═╡ e8b3920a-71ca-4dab-81af-1225516a668d
using DataFrames

# ╔═╡ c0b19231-24a1-4c77-ad66-224e00a07df4
using DelimitedFiles

# ╔═╡ 4f8b568c-3ef2-4f6c-beba-0d0ea812b193
using CairoMakie

# ╔═╡ fb1bf434-4ae3-4846-bc59-82892679b7d8
using PlutoUI

# ╔═╡ acec3106-357a-4cc3-9620-4403500248b6
using Unitful, UnitfulAstro

# ╔═╡ a45fe05e-2032-437e-ae3a-55244b761743
using PlutoUI: Slider

# ╔═╡ 55092d34-a613-487a-923b-747322777694
md"""
# Notebook to analyze Warren code
"""

# ╔═╡ b6bd7b6a-a63d-4b81-beea-6fde1230cbf2
TableOfContents()

# ╔═╡ 08919364-bb50-4311-a5ae-585c12c170fc
const data_row_predicate = !startswith(r"\s*(?:i_tm|log10)");

# ╔═╡ 5f16c34f-0564-4864-a5a5-151ca6aaa2b6
function filteredstream(filename::AbstractString; predicate = data_row_predicate)
    filtered_buffer = IOBuffer()
    infilestream = open(filename)
    for line in eachline(infilestream)
        if predicate(line)
            write(filtered_buffer, line)
            write(filtered_buffer, "\n")
        end
    end
    close(infilestream)

    seekstart(filtered_buffer)
    return filtered_buffer
end

# ╔═╡ 6db46c53-9053-4048-80b5-f2a6b4ccc03a
md"""
## Read in Fortran output data file
"""

# ╔═╡ e10c82e9-8814-44d3-8dd1-d7cb390f8c73
const data_cols = @SVector [
    :i_tm => Int, :j_ph => Int,
    :log₁₀Eᵧ_eV => Float64, :log₁₀ν_Hz => Float64,
    :log₁₀νFν_JyHz => Float64, :log₁₀Fν_Jy => Float64,
];

# ╔═╡ 593677dc-b7d5-4af1-8aa2-a45a63d4732b
dataGSdf = let
    rawdata = readdlm(filteredstream("fort.810"))
    df = DataFrame(rawdata, first.(data_cols))

    for (colname, coltype) in data_cols
        df[!, colname] = convert.(coltype, df[!, colname])
    end

    # add the exponentiated columns to the dataframe
    df.Eᵧ = exp10.(df.log₁₀Eᵧ_eV) * u"eV"
    df.ν = exp10.(df.log₁₀ν_Hz) * u"Hz"
    #df.νFν = exp10.(df.log₁₀νFν_JyHz) * u"Jy*Hz"
    df.νFν = exp10.(df.log₁₀νFν_JyHz) * 1e-26u"W/m^2"
    df.Fν = exp10.(df.log₁₀Fν_Jy) * u"Jy"

    df
end

# ╔═╡ 090e9c21-a32f-439d-a6ab-f5e1c2c66c5f


# ╔═╡ 859430ec-9177-4927-8162-f25f70df4dd1
md"""
Group up the data frame based on `i_tm`: (`**df` vs `**grouped`)
"""

# ╔═╡ ad358e63-1588-4a8a-8425-b6125034f578
dataGSgrouped = groupby(dataGSdf, :i_tm);

# ╔═╡ dd302f05-94b0-4f84-b892-b9ef4b6416fd
md"""
Also read in the `t_obs` for illustrative purposes
"""

# ╔═╡ e9e949ea-b111-4d0a-b721-38098bf6a039
t_obs = readdlm(filteredstream("fort.810", predicate = startswith("i")))[:,end]

# ╔═╡ 11253b46-7b48-4025-83a6-405454397ac5
md"""
## Plot single time index
"""

# ╔═╡ 40a1484e-4860-462a-942f-034b6d4fb9d3
md"""
Time index to plot:

`i_tm_to_plot` = $(
    @bind i_tm_to_plot Slider(axes(dataGSgrouped,1), show_value = true)
)
"""

# ╔═╡ 53d80560-b3db-4e55-ab41-e491aae70bde
md"""
Plotting at time `t_obs` = $(t_obs[i_tm_to_plot])
"""

# ╔═╡ cbc3c0db-38b2-4569-b493-6b5bab683421
lines(
    dataGSgrouped[i_tm_to_plot].log₁₀ν_Hz,
    dataGSgrouped[i_tm_to_plot].log₁₀Fν_Jy,
    axis = (;
        xlabel = "log₁₀(ν / Hz)", ylabel = "log₁₀(F_ν / Jy)",
    )
)

# ╔═╡ ce5ddcb5-ad5f-41fd-8551-9d66b8a3cced
lines(
    dataGSgrouped[i_tm_to_plot].log₁₀ν_Hz,
    dataGSgrouped[i_tm_to_plot].log₁₀νFν_JyHz,
    axis = (;
        xlabel = "log₁₀(ν / Hz)", ylabel = "log₁₀(νF_ν / Jy⋅Hz)",
    )
)

# ╔═╡ 7f5367a5-b8b8-42c0-90c3-1c06b46a625b
let df = dataGSgrouped[i_tm_to_plot]

    lines(df.ν, df.νFν, axis = (; xlabel = "ν", ylabel = "νF_ν"))
end

# ╔═╡ cbc6bff0-bd86-48a5-b9e2-6680ade0d6fb
md"""
Plot of the derivative
"""

# ╔═╡ b9b27d47-340d-4d2d-a5a4-318627ea9199
md"""
## Plot all time indices overlaid
"""

# ╔═╡ 53365cc9-9c80-42b2-bdc7-505641945ad9
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "log₁₀(F/Jy)")

    for df in dataGSgrouped[begin:2:end]
        lines!(ax, df.log₁₀ν_Hz, df.log₁₀Fν_Jy)
    end

    fig
end

# ╔═╡ 7b8d17fb-1a5e-47c2-af77-afe76e8efc03
md"""
Plot of the derivative, sort of
"""

# ╔═╡ 29fe491b-416e-40cd-91d0-d5d9d44346fb
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "log₁₀(νFν/Jy⋅Hz)")

    for df in dataGSgrouped[begin:3:end]
        lines!(ax, df.log₁₀ν_Hz, df.log₁₀νFν_JyHz)
    end

    fig
end

# ╔═╡ c2c97713-aa2e-48e8-a8ce-bbbb2e1b502e
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "ν", ylabel = "νF_ν")

    for df in dataGSgrouped[begin:3:end]
        lines!(ax, df.ν, df.νFν)
    end

    fig
end

# ╔═╡ 27eeda58-2414-4587-a5da-97f6feef6b98
midpoints(a) = (a[begin:end-1] + a[begin+1:end])/2

# ╔═╡ d1b1b7b8-97ff-49ca-a981-cc597795cd23
let
    log₁₀ν = dataGSgrouped[i_tm_to_plot].log₁₀ν_Hz
    dlog₁₀ν = diff(log₁₀ν)
    dlog₁₀Fν = diff(dataGSgrouped[i_tm_to_plot].log₁₀Fν_Jy)

    f = Figure()
    ax = Axis(f[1,1], xlabel = "log₁₀(ν / Hz)", ylabel = "dlog₁₀(F_ν / Jy)")
    lines!(ax, midpoints(log₁₀ν), dlog₁₀Fν ./ dlog₁₀ν)
    f
end

# ╔═╡ 01ae330d-56e4-4809-8f07-e6d5418dae97
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "log₁₀(ν / Hz)", ylabel = "dlog₁₀(F / Jy) / dlog₁₀(ν / Hz)")

    for df in dataGSgrouped[begin:2:end]
        log₁₀ν = df.log₁₀ν_Hz
        dlog₁₀ν = diff(log₁₀ν)
        dlog₁₀F = diff(df.log₁₀Fν_Jy)
        lines!(ax, midpoints(log₁₀ν), dlog₁₀F./dlog₁₀ν)
    end

    fig
end

# ╔═╡ 383af93b-d9cd-4c20-8a98-b80da090b5f3
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "log₁₀(ν / Hz)", ylabel = "dlog₁₀(F / Jy) − dlog₁₀(ν / Hz)")

    for df in dataGSgrouped[begin:2:end]
        log₁₀ν = df.log₁₀ν_Hz
        dlogν = diff(log₁₀ν)
        dlogF = diff(df.log₁₀Fν_Jy)
        lines!(ax, midpoints(log₁₀ν), dlogF.-dlogν)
    end

    fig
end

# ╔═╡ 4d6a5e23-9400-49e8-81be-50b36838d15d
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "ν", ylabel = "dF_ν / dν")

    for df in dataGSgrouped[begin:2:end]
        dν = diff(df.ν)
        dF = diff(df.Fν)
        scatter!(ax, midpoints(df.ν), dF ./ dν)
        #lines!(ax, midpoints(df.ν), dF .- dν)
    end

    fig
end

# ╔═╡ Cell order:
# ╟─55092d34-a613-487a-923b-747322777694
# ╠═aee7b830-cd3f-11ef-2c81-6b39b7339e04
# ╠═d96fe018-58aa-4e2f-b6d5-272f66c2ce8b
# ╠═e8b3920a-71ca-4dab-81af-1225516a668d
# ╠═c0b19231-24a1-4c77-ad66-224e00a07df4
# ╠═4f8b568c-3ef2-4f6c-beba-0d0ea812b193
# ╠═fb1bf434-4ae3-4846-bc59-82892679b7d8
# ╠═acec3106-357a-4cc3-9620-4403500248b6
# ╠═b6bd7b6a-a63d-4b81-beea-6fde1230cbf2
# ╠═5f16c34f-0564-4864-a5a5-151ca6aaa2b6
# ╠═08919364-bb50-4311-a5ae-585c12c170fc
# ╟─6db46c53-9053-4048-80b5-f2a6b4ccc03a
# ╠═593677dc-b7d5-4af1-8aa2-a45a63d4732b
# ╠═e10c82e9-8814-44d3-8dd1-d7cb390f8c73
# ╠═090e9c21-a32f-439d-a6ab-f5e1c2c66c5f
# ╟─859430ec-9177-4927-8162-f25f70df4dd1
# ╠═ad358e63-1588-4a8a-8425-b6125034f578
# ╟─dd302f05-94b0-4f84-b892-b9ef4b6416fd
# ╠═e9e949ea-b111-4d0a-b721-38098bf6a039
# ╟─11253b46-7b48-4025-83a6-405454397ac5
# ╠═a45fe05e-2032-437e-ae3a-55244b761743
# ╟─40a1484e-4860-462a-942f-034b6d4fb9d3
# ╟─53d80560-b3db-4e55-ab41-e491aae70bde
# ╠═cbc3c0db-38b2-4569-b493-6b5bab683421
# ╠═ce5ddcb5-ad5f-41fd-8551-9d66b8a3cced
# ╠═7f5367a5-b8b8-42c0-90c3-1c06b46a625b
# ╟─cbc6bff0-bd86-48a5-b9e2-6680ade0d6fb
# ╠═d1b1b7b8-97ff-49ca-a981-cc597795cd23
# ╟─b9b27d47-340d-4d2d-a5a4-318627ea9199
# ╠═53365cc9-9c80-42b2-bdc7-505641945ad9
# ╟─7b8d17fb-1a5e-47c2-af77-afe76e8efc03
# ╠═01ae330d-56e4-4809-8f07-e6d5418dae97
# ╠═383af93b-d9cd-4c20-8a98-b80da090b5f3
# ╠═4d6a5e23-9400-49e8-81be-50b36838d15d
# ╠═29fe491b-416e-40cd-91d0-d5d9d44346fb
# ╠═c2c97713-aa2e-48e8-a8ce-bbbb2e1b502e
# ╠═27eeda58-2414-4587-a5da-97f6feef6b98
