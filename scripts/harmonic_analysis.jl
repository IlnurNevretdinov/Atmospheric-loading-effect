### A Pluto.jl notebook ###
# v0.19.36


# ╔═╡ df2fa08e-bb75-11ee-2b8d-45c0a42dd9b0

using AssociatedLegendrePolynomials, SHTOOLS
using OffsetArrays, Measurements, ProgressMeter
using LinearAlgebra, Tullio, DelimitedFiles	
using Dierckx
using Rasters, NCDatasets, ArchGDAL
using CairoMakie, .Threads



# ╔═╡ 47051a23-d661-4cab-9241-72ef6982fee0
md"""
The original view of the formula to compute the harmonic coefficients $C_{nk}, S_{nk}$ of a some function:

$$\begin{equation}
	\left.\begin{array}{l}
	\bar{C}_{nk} \\
	\bar{S}_{nk}
	\end{array}\right\}=\frac{1}{4 \pi} \int \limits_\sigma f(\theta, \lambda) \bar{P}_{nk}(\sin \varphi)
        \left\{ \begin{array}{l}
	                \cos k \lambda \\
	                \sin k \lambda
	            \end{array}
        \right\} d \sigma.
\end{equation}$$

The function which to be analyzed has the view $f(x,y) = \exp(x) ⋅ \exp(y)$.

By means the package [Charm](https://github.com/blazej-bucha/charm) I defined a global grid with 
resolution $1^{\circ}×1^{\circ}$ on which the values of function were computed. Grid's 
vectors of coordinates defined as:

* latitude φ in a range __[-π/2 … π/2]__;
* longitude λ in a range __[0 … 2π]__.

Then using an approximate quadrature I calculated the harmonic coefficients.

Now I am loading all the mentioned data.
"""

load_data = Raster("../data/exp_raw/surf_pressure.nc")
# convert to hPa
# global_pressure = load_data[Ti(1), X(0.25 .. 359.75), Y(-89.75 .. 89.75)]' ./ 100
global_pressure = load_data[Ti(1)]' ./ 100

lon = Array(dims(global_pressure, X))
lat = Array(dims(global_pressure, Y))
ref_data =  Float64.(Array(global_pressure))


function pressure_heatmap(grid, data)
    x = rad2deg.(grid.first_axis)
    y = rad2deg.(grid.second_axis)
    fig = Figure(size = (1600, 800), fontsize = 25)
    ax = Axis(fig[1,1], xlabel = "λ", ylabel = "φ", )
    hmap = heatmap!(ax, x, y, data', colormap = :Spectral_11)
    # contour!(ax, x, y, data', levels=-400:100:300, color=:black, linewidth=0.85)
    Colorbar(fig[1,2], hmap, label = "hPa")
    fig
end



# # ╔═╡ 4f821328-5138-4e89-aa27-5533928aad36
# md"""
# Define a structure to store the information about global grid. It will
# have two vectors of values along each coordinate instead holding 
# the entire **2D** grid.
# """
begin
	mutable struct RegularGrid{T<:AbstractVector, U<:Tuple}
    	first_axis::T
    	second_axis::T
    	step::U
	end

	function RegularGrid(x, y)
		step = abs(x[1]-x[2]), abs(y[1] - y[2])
    	return RegularGrid(x, y, step)
	end
	
end

grid = RegularGrid(deg2rad.(lon), deg2rad.(lat))

pressure_heatmap(grid, ref_data)



# ╔═╡ 12b67bd2-521a-4110-b32e-7a115de291f8
md"""
Performing the spherical analysis/synthesis requires a lot of computations Legendre function. The package 
[AssociatedLegendrePolynomials](https://github.com/jmert/AssociatedLegendrePolynomials.jl) will help us 
with this task.

Package's 4π-normalization differs from the one commonly used in geodesy, so I must fix this by simple 
multiplication by √2 needed columns.

Function for single value $\bar{P}_{nk} (\sin φ)$ will look as follows
"""

# ╔═╡ 8d092688-6c7a-4b9b-aa27-2f11b43c2971
function correct_lgn(x, n, k)
    pnk = legendre(LegendreFourPiNorm(), n, k, x)
    if k == 0
        return pnk
    else
        if isodd(k)
            return -pnk * √2 
        else
            return pnk * √2
        end
    end
end


# ╔═╡ 6e681457-d01c-42e2-bde0-aa55f1b7bdb8
md"""
Function for all values of $\bar{P}_{nk} (\sin φ)$, it will need at the spherical synthesis stage.
"""

# ╔═╡ 799f719e-8eb2-4110-81ce-d6c6d6522b60
function correct_lgn!(out, x, N)
	legendre!(LegendreFourPiNorm(), out, N, N, x)
    @views for i in 2:N+1
		v = out[i:end,i]
        if isodd(i)
            map!(x -> √2 * x, out[i:end,i], v)
        else
            map!(x -> -√2 * x, out[i:end,i], v)
        end
    return out
    end
end



# ╔═╡ 89d9c7e5-29d5-44bc-aa0a-db1c43fc8530
md"""
And the naive implementation of numerical estimation of the initial equation for calculation of harmonic 
coefficients will look as follows
"""

# ╔═╡ f75187fb-a065-4afc-99d7-653611a921bc
function direct_analysis(grid, values, N)
	t = eltype(values)
    Cnk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    Snk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    areas = prod(grid.step) * cos.(grid.second_axis)
    @showprogress for n in 0:N
        for k in 0:n
            s₁ = zero(t)
            s₂ = zero(t)
            for i in eachindex(grid.second_axis)
				pnk = correct_lgn(sin(grid.second_axis[i]), n, k)
				area = areas[i]
				for j in eachindex(grid.first_axis)
                	sin_kλ, cos_kλ = sincos(k * grid.first_axis[j])
                	s₁ += values[i,j] * pnk * cos_kλ * area
                	s₂ += values[i,j] * pnk * sin_kλ * area
            	end
			end
        Cnk[n,k] = s₁
        Snk[n,k] = s₂
        end
    end
    return 1/4π * Cnk.parent, 1/4π * Snk.parent
end



# function optimized_analysis(grid, values, N)
# 	t = eltype(values)
#     Cnk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
#     Snk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
#     areas = prod(grid.step) * cos.(grid.second_axis)
#     weighted_values = values .* areas
#     @showprogress for n in 0:N
#         for k in 0:n
#             s₁ = zero(t)
#             s₂ = zero(t)
#             for i in eachindex(grid.second_axis)
# 				pnk = correct_lgn(sin(grid.second_axis[i]), n, k)
# 				area = areas[i]
# 				for j in eachindex(grid.first_axis)
#                 	sin_kλ, cos_kλ = sincos(k * grid.first_axis[j])
#                 	s₁ += values[i,j] * pnk * cos_kλ * area
#                 	s₂ += values[i,j] * pnk * sin_kλ * area
#             	end
# 			end
#         Cnk[n,k] = s₁
#         Snk[n,k] = s₂
#         end
#     end
#     return 1/4π * Cnk.parent, 1/4π * Snk.parent
# end



# ╔═╡ abc11262-d16c-4647-ae75-e6e68f11d9cb
md"""
Perform this analysis with $N_{\max} = 100$
"""

# ╔═╡ 48288d97-df06-4f91-aa72-5d971c9feb13
N = 360

nodes, w = SHGLQ(nothing, N)
latglq, longlq = GLQGridCoord(N) 
# interpolate reference data on GLQ grid
rs_glq = Raster(rand(X(longlq), Y(latglq)))'
itp_pressure = resample(global_pressure, to = rs_glq, size = size(rs_glq))
glq_data = Array(itp_pressure)
# show interpolated values
grid_itp = RegularGrid(deg2rad.(longlq), deg2rad.(latglq))
pressure_heatmap(grid_itp, glq_data)
# perform exact analysis
cilm = SHExpandGLQ(N, glq_data, w, nothing, nodes)
glq_cnk = cilm[1,:,:] 
glq_snk = cilm[2,:,:] 
# approximate analysis
cnk,snk = direct_analysis(grid, ref_data, N)



# ╔═╡ a145d94c-698d-4d39-adbe-fad2b7d8de3d
md"""
Now it's time for synthesis!

The initial expression looks as follows:

$$f(\theta, \lambda) = \sum_{n=0}^{N} \sum_{k=0}^n \bar{P}_{n k}(\sin \varphi)\left(\bar{C}_{nk} \cos k 
\lambda+\bar{S}_{nk} \sin k \lambda\right).$$

If change the order of summation:

$$\sum_{n=0}^{N} \sum_{k=0}^n \rightarrow \sum_{k=0}^{N} \sum_{n=k}^k,$$

we will get a new computation scheme:

$$\begin{equation}
\left.\begin{array}{c}
A_k\left(\varphi_i\right) \\
B_k\left(\varphi_i\right)
\end{array}\right\} = \sum_{n=k}^N \bar{P}_{nk}\left(\sin \varphi_i\right)\left\{\begin{array}{l}
\bar{C}_{nk} \\
\bar{S}_{nk},
\end{array}\right.
\end{equation}$$

$$\begin{equation}
f\left(\varphi_i, \lambda_j\right)=\sum^N_{k=0} A_k\left(\varphi_i\right) \cos k\lambda_j + 
B_k \left(\varphi_i\right) \sin k \lambda_j
\end{equation}$$

"""

# ╔═╡ d19cf28c-e214-4c1a-a66d-b85ff11287c4
md"""
The coefficients $A_k,B_k$ sometimes called __lumped coefficients__. Their estimation, essentially, 
leads to a series of scalar product operations, each of which is efficiently implemented in the BLAS libraries.
"""

# ╔═╡ 0b77298f-9671-433c-8a6e-2683eb34f5a4
function lumped_coefficients!(aₖ, bₖ, pnk, cnk, snk, N)
    for k in 0:N
		idx = k + 1
        @views pnk_vec = pnk[idx:end,idx] 
        @views aₖ[idx] = dot(pnk_vec, cnk[idx:end,idx])
        @views bₖ[idx] = dot(pnk_vec, snk[idx:end,idx])
    end
    return aₖ, bₖ
end

# ╔═╡ f72f11ca-8adc-4df3-9d74-92bf3adb4f50
md"""
Function to calculate array with all orders up to order __`k`__ of given `i` longitudes. 
It returns an array with size `(k+1,i)`, where each column is a vector of values `F(k⋅λ)`, 
`F` -- a function of __`cos`__ or __`sin`__.
"""

# ╔═╡ da9fb998-d90b-4143-8a83-2feafcc07522
compute_order_series(vector, order_max; f) = f.(vector' .* collect(0:order_max));


# ╔═╡ 1c5cf9a2-2445-4c44-93d8-329b5f6e0b52
md"""
A summation along $i-$th latitude performs using Einstein notation via the package
[Tullio](https://github.com/mcabbott/Tullio.jl), which takes benefit from SIMD-instructions:

$$s_i = a_k \cdot t_{k,i} + b_k \cdot q_{k,i},$$

where $t - \cos k\lambda$ and $q - \sin k\lambda$.
"""
# ╔═╡ c274c445-5d84-40e5-8485-3c217f857f76
function sum_tullio!(out, aₖ, cos_kλ, bₖ, sin_kλ)
	@tullio threads = false out[j] = aₖ[n] * cos_kλ[n,j] + bₖ[n] * sin_kλ[n,j]
end	


# ╔═╡ 5328752f-e30d-4073-83fe-3d32e2190d08
md"""
The main synthesis function for regural grid of points
"""
# ╔═╡ c7881c46-9ec2-428f-a7c6-2826d19d17dd
function synthesis!(out, grid, cnk, snk, N)
    cos_kλ = compute_order_series(grid.first_axis, N, f = cos)
    sin_kλ = compute_order_series(grid.first_axis, N, f = sin)
    pnk = zeros(N+1,N+1)
	aₖ = zeros(N+1)
    bₖ = zeros(N+1)
    @showprogress for (i,φ) in enumerate(grid.second_axis)
        correct_lgn!(pnk, sin(φ), N)
        lumped_coefficients!(aₖ, bₖ, pnk, cnk, snk, N)
        sum_tullio!(view(out, i, :), aₖ, cos_kλ, bₖ, sin_kλ)
    end
    return out
end


out_da = zeros(length(grid.second_axis), length(grid.first_axis))
out_glq = zeros(length(grid.second_axis), length(grid.first_axis))
synthesis!(out_da, grid, cnk, snk, N)
synthesis!(out_glq, grid, glq_cnk, glq_snk, N)

pressure_heatmap(grid, out_da)
pressure_heatmap(grid, out_glq)
pressure_heatmap(grid, out_da - out_glq)
# comparison with reference data
pressure_heatmap(grid, out_da - ref_data)
pressure_heatmap(grid, out_glq - ref_data)


# save results as a NetCDF file
rs_err = Raster(rand(X(lon), Y(lat)))
rs_err.data .= (out - ref_data)'
write("../data/sims/expansion_errors.nc", rs_err, force = true)




lln = readdlm("../data/exp_pro/lln_PREM.txt", skipstart = 14)
h_num = lln[1:N+1, 2]
k_num = lln[1:N+1, 4]
delta_factor = [(h_num[n+1] - (n + 1) / 2 * k_num[n+1]) / (2n + 1) for n in 0:N]


function synthesis(point, delta_factor, cnk, snk, N)
    λ,φ = point
    f = -6 / (5.5134e+3 * 6371000)
    cos_kλ = compute_order_series(λ, N, f = cos)
    sin_kλ = compute_order_series(λ, N, f = sin)

    pnk = zeros(N+1,N+1)
    correct_lgn!(pnk, sin(φ), N)
    s = sum((cnk .* cos_kλ' + snk .* sin_kλ') .* pnk, dims = 2)
    
    return f * dot(delta_factor, s)

end



p = deg2rad.([37.51604, 55.85503])

synthesis(p, delta_factor, cnk, snk, N) * 1e+9



# ╔═╡ 205be24d-df7f-4850-a60b-b044e8da8b8a
md"""
Load the data for synthesis which were generated in Charm package
"""

# ╔═╡ 365d7c7e-3498-48ce-b0f6-be348953941b
begin
	lat_mod = vec(readdlm("analysis/synthesis_data/lat_synth.txt"))
	lon_mod = vec(readdlm("analysis/synthesis_data/lon_synth.txt"))
	grid_mod = RegularGrid(lon_mod, lat_mod)
end

# ╔═╡ 0a3637db-927b-4845-abae-f9b7a17db9fe
md"""
Spherical expansion coefficients from Charm, they have to be converted in lower triangle
"""

# ╔═╡ c0609fad-00b8-4276-8227-3f220124e741
#=╠═╡
begin
	function vec_to_upper(data, n)
    	[i<=j ? data[floor(Int, j*(j-1)/2+i)] : 0. for i=1:n, j=1:n]
	end
	
	ch_coeff = readdlm("analysis/analysis_data/harmonic_coefficients.txt", skipstart=1)
	
	ch_cnk = permutedims(vec_to_upper(ch_coeff[:, 3], N+1))
	ch_snk = permutedims(vec_to_upper(ch_coeff[:, 4], N+1))
end;
  ╠═╡ =#


# ╔═╡ de555387-f00f-413b-bb9c-e912f17518be
#=╠═╡
begin
	mod_data = readdlm("analysis/synthesis_data/modelled_values.txt")
	extrema((out - mod_data) ./ mod_data * 100)
end
  ╠═╡ =#

# ╔═╡ 4668ab85-c915-4c6e-af51-613750f5181b
sh_dtm = readdlm("analysis/synthesis_data/Earth2012.RET2012.SHCto2160.dat", skipstart=1)

# ╔═╡ 0ac9a515-df58-4e02-8577-4511b0bc2858
#=╠═╡
begin
	N_dtm = 2160
	dtm_cnk = permutedims(vec_to_upper(sh_dtm[:, 3], N_dtm+1))
	dtm_snk = permutedims(vec_to_upper(sh_dtm[:, 4], N_dtm+1))
end
  ╠═╡ =#

# ╔═╡ 31b98555-7edc-4e56-9040-42a5f1607d48
begin
	step_dtm = 0.003
	lat_dtm = [40 + step_dtm - i * step_dtm for i in 1:1799]
	lon_dtm = [250 + step_dtm - i * step_dtm for i in 1:3599]
	grid_dtm = RegularGrid(deg2rad.(lon_dtm), deg2rad.(lat_dtm))
end

# ╔═╡ 394287ac-92f7-4806-9e77-2701a5056dde
begin
	out_dtm = zeros(length(lat_dtm), length(lon_dtm))
	# synthesis!(out_dtm, grid_dtm, dtm_cnk, dtm_snk, N_dtm)
end

# ╔═╡ d2f58814-5d6c-48cb-8049-3998229651e0
extrema(out_dtm)

# ╔═╡ f76937f3-f94f-464c-9245-ab6bd5db8eb0


# ╔═╡ 5f65d9a7-8542-4617-8456-a159901faa2f


# ╔═╡ eb67da2d-6d34-45f9-a5f5-e6feb8706986


# ╔═╡ 67972ca2-a952-42cb-9b5c-dea60b2e2642


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AssociatedLegendrePolynomials = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
HCubature = "19dc6840-f33b-545b-b366-655c7e3ffd49"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Tullio = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"

[compat]
AssociatedLegendrePolynomials = "~1.0.1"
DelimitedFiles = "~1.9.1"
HCubature = "~1.5.1"
Measurements = "~2.11.0"
OffsetArrays = "~1.13.0"
PlutoUI = "~0.7.55"
ProgressLogging = "~0.1.4"
Tullio = "~0.3.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "c61cc1ba91701d8b3424525c26d797766b0b617c"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssociatedLegendrePolynomials]]
git-tree-sha1 = "3204d769e06c5678b23cf928d850f2f4ad5ec8a5"
uuid = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "75bd5b6fc5089df449b5d35fa501c846c9b6549b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.12.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "e95b36755023def6ebc3d269e6483efa8b2f7f65"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.5.1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "Requires"]
git-tree-sha1 = "bdcde8ec04ca84aef5b124a17684bf3b302de00e"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.11.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "7b0e9c14c624e435076d19aea1e5cbdec2b9ca37"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.2"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.Tullio]]
deps = ["DiffRules", "LinearAlgebra", "Requires"]
git-tree-sha1 = "6d476962ba4e435d7f4101a403b1d3d72afe72f3"
uuid = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
version = "0.3.7"

    [deps.Tullio.extensions]
    TullioCUDAExt = "CUDA"
    TullioChainRulesCoreExt = "ChainRulesCore"
    TullioFillArraysExt = "FillArrays"
    TullioTrackerExt = "Tracker"

    [deps.Tullio.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═df2fa08e-bb75-11ee-2b8d-45c0a42dd9b0
# ╟─47051a23-d661-4cab-9241-72ef6982fee0
# ╠═fd10571f-1576-48b6-a26f-43c2bda80f92
# ╟─4f821328-5138-4e89-aa27-5533928aad36
# ╠═6589c417-f7a5-4812-99f5-518691acd7a7
# ╟─83a6c833-63a2-484e-a3a3-5a082e2dff7d
# ╠═e368d068-5778-4f26-bfed-3125905b6938
# ╠═b504e52c-2918-4b71-8816-a5e93ca290e7
# ╠═6f9778a2-b44a-4cbb-b6f9-6a598e713fae
# ╠═3837ddf9-5950-44c6-8a8d-ee83767e0538
# ╠═35e528c8-495c-46cb-a5a1-b46eee9a22c3
# ╟─12b67bd2-521a-4110-b32e-7a115de291f8
# ╠═8d092688-6c7a-4b9b-aa27-2f11b43c2971
# ╟─6e681457-d01c-42e2-bde0-aa55f1b7bdb8
# ╠═799f719e-8eb2-4110-81ce-d6c6d6522b60
# ╟─89d9c7e5-29d5-44bc-aa0a-db1c43fc8530
# ╠═5b79d1ef-78d0-423f-aecb-b2ef4fd1de8d
# ╠═f75187fb-a065-4afc-99d7-653611a921bc
# ╟─abc11262-d16c-4647-ae75-e6e68f11d9cb
# ╠═48288d97-df06-4f91-aa72-5d971c9feb13
# ╠═fca88aab-fb0d-426f-a2a6-f63763c787a0
# ╟─a145d94c-698d-4d39-adbe-fad2b7d8de3d
# ╟─d19cf28c-e214-4c1a-a66d-b85ff11287c4
# ╠═0b77298f-9671-433c-8a6e-2683eb34f5a4
# ╟─f72f11ca-8adc-4df3-9d74-92bf3adb4f50
# ╠═da9fb998-d90b-4143-8a83-2feafcc07522
# ╟─1c5cf9a2-2445-4c44-93d8-329b5f6e0b52
# ╠═c274c445-5d84-40e5-8485-3c217f857f76
# ╟─5328752f-e30d-4073-83fe-3d32e2190d08
# ╠═c7881c46-9ec2-428f-a7c6-2826d19d17dd
# ╟─205be24d-df7f-4850-a60b-b044e8da8b8a
# ╠═365d7c7e-3498-48ce-b0f6-be348953941b
# ╟─0a3637db-927b-4845-abae-f9b7a17db9fe
# ╠═c0609fad-00b8-4276-8227-3f220124e741
# ╠═f9b9fd74-f023-4382-acf8-91322648bbef
# ╠═de555387-f00f-413b-bb9c-e912f17518be
# ╠═4668ab85-c915-4c6e-af51-613750f5181b
# ╠═0ac9a515-df58-4e02-8577-4511b0bc2858
# ╠═31b98555-7edc-4e56-9040-42a5f1607d48
# ╠═394287ac-92f7-4806-9e77-2701a5056dde
# ╠═d2f58814-5d6c-48cb-8049-3998229651e0
# ╠═f76937f3-f94f-464c-9245-ab6bd5db8eb0
# ╠═5f65d9a7-8542-4617-8456-a159901faa2f
# ╠═eb67da2d-6d34-45f9-a5f5-e6feb8706986
# ╠═67972ca2-a952-42cb-9b5c-dea60b2e2642
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
