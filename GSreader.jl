### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

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

# ╔═╡ 40a1484e-4860-462a-942f-034b6d4fb9d3
i_tm_to_plot = 10;

# ╔═╡ e10c82e9-8814-44d3-8dd1-d7cb390f8c73
const data_cols = @SVector [
	:i_tm => Int, :j_ph => Int,
	:log10eV => Float64, :log10Hz => Float64,
	:log10νFν => Float64, :log10Jy => Float64,
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

# ╔═╡ cbc3c0db-38b2-4569-b493-6b5bab683421
lines(dataGSgdf[i_tm_to_plot].log10Hz, dataGSgdf[i_tm_to_plot].log10Jy)

# ╔═╡ 53365cc9-9c80-42b2-bdc7-505641945ad9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log10(ν/Hz)", ylabel = "log10(F/Jy)")

	for df in dataGSgdf[begin:2:end]
		lines!(ax, df.log10Hz, df.log10Jy)
	end

	fig
end

# ╔═╡ 29fe491b-416e-40cd-91d0-d5d9d44346fb
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log10(ν/Hz)", ylabel = "log10(νFν/Jy⋅Hz)")

	for df in dataGSgdf[begin:3:end]
		lines!(ax, df.log10Hz, df.log10νFν)
	end

	fig
end

# ╔═╡ Cell order:
# ╠═aee7b830-cd3f-11ef-2c81-6b39b7339e04
# ╠═d96fe018-58aa-4e2f-b6d5-272f66c2ce8b
# ╠═e8b3920a-71ca-4dab-81af-1225516a668d
# ╠═c0b19231-24a1-4c77-ad66-224e00a07df4
# ╠═4f8b568c-3ef2-4f6c-beba-0d0ea812b193
# ╠═5f16c34f-0564-4864-a5a5-151ca6aaa2b6
# ╠═08919364-bb50-4311-a5ae-585c12c170fc
# ╠═593677dc-b7d5-4af1-8aa2-a45a63d4732b
# ╠═ad358e63-1588-4a8a-8425-b6125034f578
# ╠═40a1484e-4860-462a-942f-034b6d4fb9d3
# ╠═cbc3c0db-38b2-4569-b493-6b5bab683421
# ╠═53365cc9-9c80-42b2-bdc7-505641945ad9
# ╠═29fe491b-416e-40cd-91d0-d5d9d44346fb
# ╠═e10c82e9-8814-44d3-8dd1-d7cb390f8c73
