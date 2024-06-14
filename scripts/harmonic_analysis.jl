using SHTOOLS
using DelimitedFiles	
using Rasters, NCDatasets, ArchGDAL
using CairoMakie, .Threads



# auxiliary functions
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

mutable struct RegularGrid{T<:AbstractVector, U<:Tuple}
    first_axis::T
    second_axis::T
    step::U
end

function RegularGrid(x, y)
    step = abs(x[1]-x[2]), abs(y[1] - y[2])
    return RegularGrid(x, y, step)
end



include("../src/spharm_utils.jl")

# loading data about atmospheric surface pressure
load_data = Raster("../data/exp_raw/surf_pressure.nc")
# global_pressure = load_data[Ti(1), X(0.25 .. 359.75), Y(-89.75 .. 89.75)]' ./ 100
global_pressure = load_data[Ti(1)]' ./ 100 # convert to hPa

lon = Array(dims(global_pressure, X))
lat = Array(dims(global_pressure, Y))
grid = RegularGrid(deg2rad.(lon), deg2rad.(lat))
ref_data =  Float64.(Array(global_pressure))
# image of data
pressure_heatmap(grid, ref_data)


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
cnk,snk = sphharm.direct_analysis(grid, ref_data, N)

# synthesis with harmonic coefficients from direct analysis
out_da = zeros(length(grid.second_axis), length(grid.first_axis))
sphharm.synthesis!(out_da, grid, cnk, snk, N)
pressure_heatmap(grid, out_da)

# the same but with GLQ-based solution
out_glq = zeros(length(grid.second_axis), length(grid.first_axis))
sphharm.synthesis!(out_glq, grid, glq_cnk, glq_snk, N)
pressure_heatmap(grid, out_glq)

# difference
pressure_heatmap(grid, out_da - out_glq)

# comparison both solution with reference data
pressure_heatmap(grid, out_da - ref_data)
pressure_heatmap(grid, out_glq - ref_data)


# save results as a NetCDF file
# rs_err = Raster(rand(X(lon), Y(lat)))
# rs_err.data .= (out - ref_data)'
# write("../data/sims/expansion_errors.nc", rs_err, force = true)


# Love numbers from LoadDef
lln = readdlm("../data/exp_pro/lln_PREM.txt", skipstart = 14)
h_num = lln[1:N+1, 2]
k_num = lln[1:N+1, 4]
delta_factor = [(h_num[n+1] - (n + 1) / 2 * k_num[n+1]) / (2n + 1) for n in 0:N]

# Moscow's coordinates
p = deg2rad.([37.51604, 55.85503])
sphharm.synthesis(p, delta_factor, cnk, snk, N) * 1e+9


