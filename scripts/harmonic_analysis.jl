using SHTOOLS
using DelimitedFiles, Statistics	
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


function initial_pressure(P, g₀, z, R, T₀)
    out = P / exp((-g₀ * z) / (R * T₀))
    return out
end



include("../src/spharm_utils.jl")

# loading data about atmospheric surface pressure
i = 1
hourly_pressure = Raster("../data/exp_raw/surf_pressure.nc")
global_pressure = hourly_pressure[Ti(i)]' ./ 100 # convert to hPa
ref_data =  Float64.(Array(global_pressure))
# pressure_heatmap(grid, ref_data)

# make a regular grid for global data
lon = Array(dims(global_pressure, X))
lat = Array(dims(global_pressure, Y))
grid = RegularGrid(deg2rad.(lon), deg2rad.(lat))

# geoid heights from EGM96
geoid_data = Raster("../data/exp_pro/geoid_egm96.nc")
geoid_heigts = Array(reverse(geoid_data', dims = 1)[:, 1:end-1])

# temperature from ECMWF
hourly_temperature = Raster("../data/exp_raw/temperature.nc")
global_temperature = hourly_temperature[Ti(i)]' 
temp_data = Float64.(Array(global_temperature))

# reduce the pressure values to the same surface
reduced_pressure = initial_pressure.(ref_data, 9.8, geoid_heigts, 287.04, temp_data)
sphere_pressure = deepcopy(global_pressure)
sphere_pressure.data .= reduced_pressure

# define expansion degree
N = 180
nodes, w = SHGLQ(nothing, N)
latglq, longlq = GLQGridCoord(N) 

# interpolate reduced data on GLQ grid
rs_glq = Raster(rand(X(longlq), Y(latglq)))'
itp_pressure = resample(sphere_pressure, to = rs_glq, size = size(rs_glq))
glq_data = Array(itp_pressure)
# show interpolated values
grid_itp = RegularGrid(deg2rad.(longlq), deg2rad.(latglq))
# pressure_heatmap(grid_itp, glq_data)
# perform exact analysis
cilm = SHExpandGLQ(N, glq_data, w, nothing, nodes)
glq_cnk = cilm[1,:,:] 
glq_snk = cilm[2,:,:] 

# approximate analysis
cnk,snk = sphharm.direct_analysis(grid, sphere_pressure, N)

q,w = sphharm.optimized_analysis(grid, sphere_pressure, N)
# @profview sphharm.optimized_analysis(grid, ref_data, 60)


# synthesis with harmonic coefficients from direct analysis
out_da = zeros(length(grid.second_axis), length(grid.first_axis))
sphharm.synthesis!(out_da, grid, cnk, snk, N) 
pressure_heatmap(grid, out_da) 

# comparison with reference data
error_da = out_da - reduced_pressure
pressure_heatmap(grid, error_da)  
extrema(error_da), mean(abs, error_da) 

# the same but for GLQ-based solution
out_glq = zeros(length(grid.second_axis), length(grid.first_axis))
sphharm.synthesis!(out_glq, grid, glq_cnk, glq_snk, N)
pressure_heatmap(grid, out_glq) 

# comparison with reference data
error_glq = out_glq - zero_surface
pressure_heatmap(grid, error_glq) 
extrema(error_glq), mean(abs, error_glq)

# difference between solutions
pressure_heatmap(grid, out_da - out_glq)


# save results as a NetCDF file
# rs_err = Raster(rand(X(lon), Y(lat)))
# rs_err.data .= (out - ref_data)'
# write("../data/sims/expansion_errors.nc", rs_err, force = true)


# Love numbers from LoadDef
lln = readdlm("../data/exp_pro/lln_PREM.txt", skipstart = 14)
h_num = lln[1:N+1, 2]
k_num = lln[1:N+1, 4] ./ lln[1:N+1,1]
k_num[1] = 0.0
delta_factor = [(h_num[n+1] - (n + 1) / 2 * k_num[n+1]) / (2n + 1) for n in 0:N]

# Moscow's coordinates
p = deg2rad.([37.51604, 55.85503])
sphharm.synthesis(p, delta_factor, glq_cnk, glq_snk, N)


