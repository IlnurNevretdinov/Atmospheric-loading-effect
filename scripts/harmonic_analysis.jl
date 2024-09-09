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


# Love numbers from LoadDef
lln = readdlm("../data/exp_pro/lln_PREM.txt", skipstart = 14)
h_num = lln[1:N+1, 2]
k_num = lln[1:N+1, 4] ./ collect(0:N)
k_num[1] = 0.0
delta_factor = [(h_num[n+1] - (n + 1) / 2 * k_num[n+1]) / (2n + 1) for n in 0:N]
# t = sphharm.f * 10^11 * delta_factor 

# Moscow's coordinates
p = deg2rad.([37.51604, 55.85503])

# define expansion degree
N = 180

# loading data about atmospheric surface pressure
hourly_pressure = Raster("../data/exp_raw/surf_pressure.nc")
# temperature from ECMWF
hourly_temperature = Raster("../data/exp_raw/temperature.nc")


# make a regular grid for global data
lon = Array(dims(hourly_pressure, X))
lat = Array(dims(hourly_pressure, Y))
grid = RegularGrid(deg2rad.(lon), deg2rad.(lat))

# geoid heights from EGM96
geoid_data = Raster("../data/exp_pro/geoid_egm96.nc")
geoid_heigts = Array(reverse(geoid_data', dims = 1)[:, 1:end-1])


begin
    i = 1
    global_pressure = hourly_pressure[Ti(i)]' ./ 100 # convert to hPa
    ref_data = Float64.(Array(global_pressure))

    global_temperature = hourly_temperature[Ti(i)]'  
    temp_data = Float64.(Array(global_temperature))
    
    # reduce all pressure values to the same surface
    reduced_pressure = initial_pressure.(ref_data, 9.8, geoid_heigts, 287.04, temp_data)
    sphere_pressure = deepcopy(global_pressure)
    sphere_pressure.data .= reduced_pressure .- mean(reduced_pressure)
end

# GLQ solution
begin
    nodes, w = SHGLQ(nothing, N)
    latglq, longlq = GLQGridCoord(N) 
    # interpolate reduced atmospheric data on GLQ grid
    rs_glq = Raster(rand(X(longlq), Y(latglq)))'
    itp_pressure = resample(sphere_pressure, to = rs_glq, size = size(rs_glq))
    glq_data = Array(itp_pressure)
    # grid_itp = RegularGrid(deg2rad.(longlq), deg2rad.(latglq))
    # pressure_heatmap(grid_itp, glq_data)
    # perform exact analysis
    cilm = SHExpandGLQ(N, glq_data, w, nothing, nodes)
    glq_cnk = cilm[1,:,:] 
    glq_snk = cilm[2,:,:] 
end

# synthesis with harmonic coefficients from GLQ-based solution
begin
    # comparison on GLQ grid
    gridglq = MakeGridGLQC(cilm, N, nothing, nodes)
    extrema(gridglq - glq_data), mean(abs, gridglq - glq_data)
    # comparison with reference data on source grid
    out_glq = zeros(length(grid.second_axis), length(grid.first_axis))
    sphharm.synthesis!(out_glq, grid, glq_cnk, glq_snk, N)
    error_glq = out_glq - sphere_pressure.data
    extrema(error_glq), mean(abs, error_glq)
    # pressure_heatmap(grid, error_glq) 
    # pressure_heatmap(grid, out_glq) 
end

sphharm.synthesis(p, delta_factor, glq_cnk, glq_snk, N) * 10^10


# approximate analysis
cnk,snk = sphharm.direct_analysis(grid, sphere_pressure.data, N)

q,w = sphharm.optimized_analysis(grid, sphere_pressure.data, N)
# @profview sphharm.optimized_analysis(grid, ref_data, 60)

# synthesis with harmonic coefficients from direct analysis
begin
    out_da = zeros(length(grid.second_axis), length(grid.first_axis))
    sphharm.synthesis!(out_da, grid, q, w, N) 
    # comparison with reference data
    error_da = out_da - sphere_pressure.data
    extrema(error_da), mean(abs, error_da) 

    # pressure_heatmap(grid, error_da)  
    # pressure_heatmap(grid, out_da) 
end

# difference between solutions
pressure_heatmap(grid, out_da - out_glq)

# save results as a NetCDF file
# rs_err = Raster(rand(X(lon), Y(lat)))
# rs_err.data .= (out - ref_data)'
# write("../data/sims/expansion_errors.nc", rs_err, force = true)





