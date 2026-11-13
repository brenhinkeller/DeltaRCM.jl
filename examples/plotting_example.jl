# Plotting a DeltaRCM model with Makie, and animating a delta as it grows.
#
#   julia --project=examples examples/demo.jl
#
# (first run: `julia --project=examples -e 'using Pkg; Pkg.instantiate()'`)
#
# Writes to `examples/output/`:
#   delta_final.png   one-line snapshot of the bed elevation from `plot_field`
#   delta_panels.png  a hand-built multi-panel figure of the final state
#   delta_growth.gif  an animation of the delta growing (bed, flow speed, surface sand fraction)
#
# Loading a Makie backend is all it takes to enable plotting: the plotting methods live in
# a package extension that activates when `Makie` is loaded.

using DeltaRCM
using CairoMakie

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

## --- a model ------------------------------------------------------------------------------
# A 3 km x 6 km basin on a 50 m grid (60 x 120 cells) grows a recognizable delta in a few
# hundred steps. Drop the keywords to use the pyDeltaRCM defaults (5 km x 10 km).
model = DeltaModel(; Length = 3000, Width = 6000, dx = 50, Np_water = 1000, Np_sed = 1000, seed = 42)

nsteps = 1500         # time steps to run
steps_per_frame = 5   # time steps between animation frames

days(m) = m.time / 86400   # model time in days

## --- Axes helpers ---------------------------------------------------------------------------
# Model fields are (L, W) arrays whose first index runs downstream. A heatmap wants x = cross-stream,
# y = downstream, so we transpose, and reverse the y axis to put the inlet at the top.
L, W = size(model)
xs = (0:W-1) .* model.dx ./ 1000     # cross-stream, km
ys = (0:L-1) .* model.dx ./ 1000     # downstream, km
flip(A) = Float32.(permutedims(A))

## --- the animation ----------------------------------------------------------------------------
# Observables hold the plotted data; updating them redraws the figure. `record` captures each
# redraw as a frame.
eta   = Observable(flip(model.eta))
speed = Observable(flip(model.uw))
sandf = Observable(flip(model.sand_frac))
title = Observable("t = 0.0 days")

fig = Figure(size = (1300, 520))
Label(fig[0, 1:6], title; fontsize = 20)

panel(col, obs, name; colormap, colorrange) = begin
    ax = Axis(fig[1, col]; title = name, xlabel = "cross-stream (km)", ylabel = col == 1 ? "downstream (km)" : "",
              yreversed = true, aspect = DataAspect())
    col > 1 && hideydecorations!(ax; grid = false)
    hm = heatmap!(ax, xs, ys, obs; colormap, colorrange)
    Colorbar(fig[1, col+1][1, 1], hm; width = 12)
    return ax
end

# Fixed colour ranges, so frames are comparable: the basin starts at -hb and land sits near 0.
panel(1, eta, "bed elevation (m)"; colormap = :topo, colorrange = (-model.hb - 0.5, 1.0))
panel(3, speed, "flow speed (m/s)"; colormap = :viridis, colorrange = (0, model.u0))
panel(5, sandf, "surface sand fraction"; colormap = :YlOrBr, colorrange = (0, 1))
colgap!(fig.layout, 8)

function update!(model)
    eta[] = flip(model.eta)
    speed[] = flip(model.uw)
    sandf[] = flip(model.sand_frac)
    title[] = "t = $(round(days(model); digits = 1)) days   (step $(model.iter))"
end

gif = joinpath(outdir, "delta_growth.gif")
@info "Running $nsteps steps and recording $gif"
nframes = nsteps ÷ steps_per_frame
elapsed = @elapsed record(fig, gif, 0:nframes; framerate = 15, px_per_unit = 1) do frame
    frame > 0 && run!(model, steps_per_frame)
    update!(model)
end
@info "Done in $(round(elapsed; digits = 1)) s (model steps and rendering together)"

## --- snapshots of the final state ------------------------------------------------------------------
# `plot_field` is the one-liner: a heatmap of any model field, labelled with the model time.
save(joinpath(outdir, "delta_final.png"), plot_field(model, :eta))

# A custom figure from the same arrays: the final bed together with the water depth and the
# net deposition since the start (erosion, mainly in the channels, plots as the palest colour).
fig2 = Figure(size = (1300, 520))
for (col, (name, data, cmap, range)) in enumerate((
        ("bed elevation (m)", model.eta, :topo, (-model.hb - 0.5, 1.0)),
        ("water depth (m)", model.depth, :deep, (0, model.hb)),
        ("deposition since start (m)", model.eta .- model.eta_init, :amp, (0, model.hb))))
    ax = Axis(fig2[1, 2col-1]; title = name, xlabel = "cross-stream (km)", ylabel = col == 1 ? "downstream (km)" : "",
              yreversed = true, aspect = DataAspect())
    col > 1 && hideydecorations!(ax; grid = false)
    hm = heatmap!(ax, xs, ys, flip(data); colormap = cmap, colorrange = range)
    Colorbar(fig2[1, 2col], hm; width = 12)
end
Label(fig2[0, 1:6], "after $(round(days(model); digits = 1)) days"; fontsize = 20)
save(joinpath(outdir, "delta_panels.png"), fig2)

@info "Wrote" outdir readdir(outdir)

## --- End of File