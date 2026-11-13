# Saving model output to NetCDF, and reading it back.
#
#   julia --project=examples examples/netcdf_output.jl
#
# (first run: `julia --project=examples -e 'using Pkg; Pkg.instantiate()'`)
#
# Writes to `examples/output/netcdf_demo/`:
#   pyDeltaRCM_output.nc   the model output, in pyDeltaRCM's file layout
#   netcdf_timeseries.pdf  quantities computed from the file: land area and bed volume vs time
#
# Loading NCDatasets is all it takes to enable NetCDF output (it is implemented in a package
# extension). The `save_*` options and `save_dt` are the same as in pyDeltaRCM, so the same
# YAML configuration files work.

using DeltaRCM
using NCDatasets
using CairoMakie     # only used for the figure at the end

outdir = joinpath(@__DIR__, "output", "netcdf_demo")
rm(outdir; force = true, recursive = true)

## --- run a model that saves grids ------------------------------------------------------------
# `save_dt` is in model seconds. Grids are written at t = 0 and then whenever at least `save_dt`
# has elapsed since the last save. With a time step of 25 000 s, 1.5 days (129 600 s) means every
# 6th step, so the 300 steps below end exactly on a saved record.
model = DeltaModel(;
    Length = 3000, Width = 6000, dx = 50, Np_water = 1000, Np_sed = 1000, seed = 42,
    out_dir = outdir,
    save_dt = 1.5 * 86400,
    save_eta_grids = true,         # bed elevation
    save_depth_grids = true,       # water depth
    save_velocity_grids = true,    # flow speed
    save_sandfrac_grids = true,    # surface sand fraction
    save_metadata = true,          # grid geometry and parameters (implied by any save_*_grids)
)
run!(model; timesteps = 300)
finalize!(model)                   # closes the file

## --- what is in the file --------------------------------------------------------------------------
file = joinpath(outdir, "pyDeltaRCM_output.nc")
ds = NCDataset(file)
println(ds)                        # dimensions, variables and attributes at a glance

# Grids are stored as (seconds, x, y) as seen from NetCDF/Python/xarray, with x downstream and y
# cross-stream. NCDatasets lists dimensions in Julia order, so a grid reads back as an
# (y, x, time) array. Transposing a snapshot gives the layout of `model.eta`: (L, W).
@show size(ds["eta"]) dimnames(ds["eta"])
final_eta = permutedims(ds["eta"][:, :, end])
println("last record matches the final model state (the file stores Float32): ", maximum(abs.(final_eta .- model.eta)) < 1e-4)

## --- reading metadata and computing things from the output ----------------------------------------
t = ds["seconds"][:] ./ 86400                  # model time of each record, days
aux = ds.group["auxdata"]                      # fixed metadata (called `meta` with `legacy_netcdf = true`)
dx = Float64(aux["dx"][])
hb = Float64(aux["hb"][])
cell_type = permutedims(aux["cell_type"][:, :])   # -2 wall, -1 edge, 0 ocean, 1 channel
println("inlet: L0 = ", aux["L0"][], " cells long, N0 = ", aux["N0"][], " cells wide; centre column CTR = ",
        aux["CTR"][], " (0-based, as in pyDeltaRCM)")

nt = length(t)
land_km2 = zeros(nt)
volume_Mm3 = zeros(nt)
eta0 = permutedims(ds["eta"][:, :, 1])
for k in 1:nt
    eta = permutedims(ds["eta"][:, :, k])
    land_km2[k] = count((eta .> -hb / 2) .& (cell_type .== 0)) * dx^2 / 1e6
    volume_Mm3[k] = sum((eta .- eta0)[cell_type .!= -2]) * dx^2 / 1e6
end
println(nt, " records; land area ", round(land_km2[end]; digits = 2), " km², deposited volume ",
        round(volume_Mm3[end]; digits = 2), " Mm³ after ", round(t[end]; digits = 1), " days")

## --- plot it ----------------------------------------------------------------------------------------
fig = Figure(size = (1100, 420))
ax1 = Axis(fig[1, 1]; xlabel = "model time (days)", ylabel = "land area (km²)", title = "land area")
lines!(ax1, t, land_km2; linewidth = 3)
ax2 = Axis(fig[1, 2]; xlabel = "model time (days)", ylabel = "deposited volume (Mm³)", title = "bed volume gained")
lines!(ax2, t, volume_Mm3; linewidth = 3)
ax3 = Axis(fig[1, 3]; xlabel = "cross-stream (km)", ylabel = "downstream (km)", title = "final bed elevation (m)",
           yreversed = true, aspect = DataAspect())
L, W = size(model)
hm = heatmap!(ax3, (0:W-1) .* dx ./ 1000, (0:L-1) .* dx ./ 1000, permutedims(final_eta); colormap = :topo,
              colorrange = (-hb - 0.5, 1.0))
Colorbar(fig[1, 4], hm; width = 12)
save(joinpath(outdir, "netcdf_timeseries.pdf"), fig)
display(fig)

close(ds)
@info "Wrote" outdir readdir(outdir)

## To read the same file from Python:
#   import xarray as xr
#   ds = xr.open_dataset("pyDeltaRCM_output.nc")                       # eta has dims (seconds, x, y)
#   aux = xr.open_dataset("pyDeltaRCM_output.nc", group="auxdata")     # L0, N0, CTR, dx, h0, hb, cell_type, ...
