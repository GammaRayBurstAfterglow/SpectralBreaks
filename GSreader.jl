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

# ╔═╡ e10c82e9-8814-44d3-8dd1-d7cb390f8c73
const data_cols = @SVector [
	:i_tm => Int, :j_ph => Int,
	:log₁₀eV => Float64, :log₁₀Hz => Float64,
	:log₁₀νFν => Float64, :log₁₀Jy => Float64,
];

# ╔═╡ 593677dc-b7d5-4af1-8aa2-a45a63d4732b
dataGSdf = let
	rawdata = readdlm(filteredstream("fort.810"))
	df = DataFrame(rawdata, first.(data_cols))

	for (colname, coltype) in data_cols
		df[!, colname] = convert.(coltype, df[!, colname])
	end

	df
end

# ╔═╡ ad358e63-1588-4a8a-8425-b6125034f578
dataGSgdf = groupby(dataGSdf, :i_tm);

# ╔═╡ 40a1484e-4860-462a-942f-034b6d4fb9d3
@bind i_tm_to_plot NumberField(axes(dataGSgdf,1))

# ╔═╡ cbc3c0db-38b2-4569-b493-6b5bab683421
lines(dataGSgdf[i_tm_to_plot].log₁₀Hz, dataGSgdf[i_tm_to_plot].log₁₀Jy)

# ╔═╡ 53365cc9-9c80-42b2-bdc7-505641945ad9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "log₁₀(F/Jy)")

	for df in dataGSgdf[begin:2:end]
		lines!(ax, df.log₁₀Hz, df.log₁₀Jy)
	end

	fig
end

# ╔═╡ 29fe491b-416e-40cd-91d0-d5d9d44346fb
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "log₁₀(νFν/Jy⋅Hz)")

	for df in dataGSgdf[begin:3:end]
		lines!(ax, df.log₁₀Hz, df.log₁₀νFν)
	end

	fig
end

# ╔═╡ 27eeda58-2414-4587-a5da-97f6feef6b98
midpoints(a) = (a[begin:end-1] + a[begin+1:end])/2

# ╔═╡ d1b1b7b8-97ff-49ca-a981-cc597795cd23
let
	ν = dataGSgdf[i_tm_to_plot].log₁₀Hz
	dlogν = diff(ν)
	dlogF = diff(dataGSgdf[i_tm_to_plot].log₁₀Jy)
	lines(midpoints(ν), dlogF./dlogν)
end

# ╔═╡ 01ae330d-56e4-4809-8f07-e6d5418dae97
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "dlog₁₀(F/Jy) / dlog₁₀(ν/Hz)")

	for df in dataGSgdf[begin:2:end]
		ν = df.log₁₀Hz
		dlogν = diff(ν)
		dlogF = diff(df.log₁₀Jy)
		lines!(ax, midpoints(ν), dlogF./dlogν)
	end

	fig
end

# ╔═╡ 383af93b-d9cd-4c20-8a98-b80da090b5f3
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log₁₀(ν/Hz)", ylabel = "dlog₁₀(F/Jy) − dlog₁₀(ν/Hz)")

	for df in dataGSgdf[begin:2:end]
		ν = df.log₁₀Hz
		dlogν = diff(ν)
		dlogF = diff(df.log₁₀Jy)
		lines!(ax, midpoints(ν), dlogF.-dlogν)
	end

	fig
end

# ╔═╡ Cell order:
# ╠═aee7b830-cd3f-11ef-2c81-6b39b7339e04
# ╠═d96fe018-58aa-4e2f-b6d5-272f66c2ce8b
# ╠═e8b3920a-71ca-4dab-81af-1225516a668d
# ╠═c0b19231-24a1-4c77-ad66-224e00a07df4
# ╠═4f8b568c-3ef2-4f6c-beba-0d0ea812b193
# ╠═fb1bf434-4ae3-4846-bc59-82892679b7d8
# ╠═5f16c34f-0564-4864-a5a5-151ca6aaa2b6
# ╠═08919364-bb50-4311-a5ae-585c12c170fc
# ╠═593677dc-b7d5-4af1-8aa2-a45a63d4732b
# ╠═ad358e63-1588-4a8a-8425-b6125034f578
# ╠═40a1484e-4860-462a-942f-034b6d4fb9d3
# ╠═cbc3c0db-38b2-4569-b493-6b5bab683421
# ╠═d1b1b7b8-97ff-49ca-a981-cc597795cd23
# ╠═53365cc9-9c80-42b2-bdc7-505641945ad9
# ╠═01ae330d-56e4-4809-8f07-e6d5418dae97
# ╠═383af93b-d9cd-4c20-8a98-b80da090b5f3
# ╠═29fe491b-416e-40cd-91d0-d5d9d44346fb
# ╠═e10c82e9-8814-44d3-8dd1-d7cb390f8c73
# ╠═27eeda58-2414-4587-a5da-97f6feef6b98
