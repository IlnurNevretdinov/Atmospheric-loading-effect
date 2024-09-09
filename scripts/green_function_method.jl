using DelimitedFiles, Statistics	
using Rasters, NCDatasets, ArchGDAL
using CairoMakie, .Threads
using Interpolations
import GeodeticIntegrals as lib


# auxiliary functions
function create_grid(x_vec, y_vec)
    n = length(y_vec)
    k = length(x_vec)
    x = x_vec' .* ones(n)
    y = reverse(y_vec) .* ones(n,k)
    return x,y
end

function grid_from_raster(data::Raster; radians::Bool)
    y_vec = sort(Array(dims(data,Y)))
    x_vec = sort(Array(dims(data,X)))
    x,y = create_grid(x_vec, y_vec)
    step = abs(x_vec[1] - x_vec[2]), abs(y_vec[1] - y_vec[2])
    if radians
        out = RegularGrid(deg2rad.(x), deg2rad.(y), deg2rad.(step))
    else
        out = RegularGrid(x, y, step)
    end
    return out
end


include("../src/spharm_utils.jl")


# loading data about atmospheric surface pressure
i = 1
hourly_pressure = Raster("../data/exp_raw/surf_pressure.nc")
global_pressure = hourly_pressure[Ti(i)]' ./ 100 # convert to hPa
data = Float64.(Array(global_pressure))


# make a regular grid for global data
lon = deg2rad.(Array(dims(global_pressure, X)))
lat = deg2rad.(Array(dims(global_pressure, Y)))
lon_glob,lat_glob = create_grid(lon,reverse(lat))
grid = lib.RegularGrid(lon_glob, lat_glob, deg2rad.((0.25,0.25)))

# Love numbers
N = 10_000
lln = readdlm("../data/exp_pro/lln_PREM.txt", skipstart = 14)
h_num = lln[1:N+1, 2]
k_num = lln[1:N+1, 4] ./ lln[1:N+1,1]
k_num[1] = 0.0
lovenumber = lib.LoveNumber(h = h_num, l = rand(10), k = k_num)

# Moscow's coordinates
p = permutedims(deg2rad.([37.51604 55.85503]))

# global integration
ψ₀ = 1
kernel = lib.GreenVertical(truncation_radius = ψ₀, LN = lovenumber, quadrature_order = 3, degree = N)
atl_global = lib.integration(data, grid, kernel)
# interpolate on given point
itp = interpolate((reverse(lat), lon), reverse(atl_global, dims = 1), Gridded(Linear()))
itp(p[2], p[1])


# local integration
atl_local = lib.integration(p, data, grid, kernel, number_neighbors = 5) 


out_local = []
for k in 1:3:24
    gp = hourly_pressure[Ti(k)]' ./ 100 # convert to hPa
    d = Float64.(Array(gp))
    atl_local = lib.integration(p, d, grid, kernel, number_neighbors = 10) 
    push!(out_local, atl_local[1])
end




begin
    x = rad2deg.(grid.first_axis[1,:])
    y = rad2deg.(grid.second_axis[:,1])
    fig = Figure(size = (1600, 800), fontsize = 25)
    ax = Axis(fig[1,1], xlabel = "λ", ylabel = "φ", )
    hmap = heatmap!(ax, x, y, atl' * 10^18, colormap = :Spectral_11)
    # contour!(ax, x, y, data', levels=-400:100:300, color=:black, linewidth=0.85)
    Colorbar(fig[1,2], hmap)
    fig
end