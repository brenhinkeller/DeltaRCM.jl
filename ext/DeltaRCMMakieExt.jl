module DeltaRCMMakieExt

using DeltaRCM
using Makie
import DeltaRCM: plot_field, save_figure

const COLORMAPS = Dict(:eta => :topo, :stage => :viridis, :depth => :deep, :qw => :viridis, :uw => :viridis,
                       :qs => :viridis, :sand_frac => :YlOrBr)

function DeltaRCM.plot_field(m, field::Symbol = :eta; colormap = get(COLORMAPS, field, :viridis),
                             title = String(field), kwargs...)
    data = getproperty(m, field)
    L, W = size(m)
    fig = Figure(size = (800, max(300, round(Int, 800 * L / W) + 120)))
    ax = Axis(fig[1, 1]; title = "$title  (t = $(round(m.time / 86400; digits = 2)) d)", xlabel = "cross-stream (m)",
              ylabel = "downstream (m)", yreversed = true, aspect = DataAspect())
    hm = heatmap!(ax, (0:W-1) .* m.dx, (0:L-1) .* m.dx, Float32.(permutedims(data)); colormap, kwargs...)
    Colorbar(fig[1, 2], hm)
    return fig
end

function DeltaRCM.save_figure(m, field::Symbol, idx::Integer)
    p = m.params
    name = p.save_figs_sequential ? "$(field)_$(lpad(idx, 5, '0')).png" : "$(field).png"
    save(joinpath(p.out_dir, name), plot_field(m, field))
    return nothing
end

end
