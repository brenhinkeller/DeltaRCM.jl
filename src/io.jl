# Output and checkpointing.
#
# The core package decides *what* is saved and *when* (the pyDeltaRCM `save_*` options,
# `save_dt`, `checkpoint_dt`); the file formats are implemented by package extensions:
#
#   NCDatasets -> `open_netcdf!`, `write_netcdf!`, `close_netcdf!`
#   JLD2       -> `save_checkpoint`, `load_checkpoint!`

"""
    NetCDFOutput

Handle to the NetCDF output file (`pyDeltaRCM_output.nc`, in pyDeltaRCM's "sandsuet 1.0.0"
layout). Created by the model when any `save_*_grids`/`save_metadata` option is set;
requires `using NCDatasets`.
"""
mutable struct NetCDFOutput
    path::String
    dataset::Any
end

"""Gridded variables that can be written, as `name => (model field, units, long_name)` in the sandsuet schema."""
const GRID_OUTPUTS = (
    eta = (:eta, "meters", "channel_bottom__elevation"),
    stage = (:stage, "meters", "channel_water_surface__elevation"),
    depth = (:depth, "meters", "channel_water__thickness"),
    discharge = (:qw, "cubic meters per second", "channel_water_flowing__volume_rate"),
    velocity = (:uw, "meters per second", "channel_water_flowing__speed"),
    sedflux = (:qs, "cubic meters per second", "channel_water_sediment_flowing__volume_rate"),
    sandfrac = (:sand_frac, "fraction", "channel_bottom_sediment_sand__volume_fraction"),
    discharge_x = (:qx, "cubic meters per second", "channel_water_flowing__x_component_of_volume_rate"),
    discharge_y = (:qy, "cubic meters per second", "channel_water_flowing__y_component_of_volume_rate"),
    velocity_x = (:ux, "meters per second", "channel_water_flowing__x_component_of_speed"),
    velocity_y = (:uy, "meters per second", "channel_water_flowing__y_component_of_speed"),
)

"""
    saved_grids(params) -> Vector{Symbol}

Names (keys of `GRID_OUTPUTS`) of the gridded variables selected by the `save_*_grids` options.
"""
function saved_grids(p::DeltaParameters)
    sel = Symbol[]
    p.save_eta_grids && push!(sel, :eta)
    p.save_stage_grids && push!(sel, :stage)
    p.save_depth_grids && push!(sel, :depth)
    p.save_discharge_grids && push!(sel, :discharge)
    p.save_velocity_grids && push!(sel, :velocity)
    p.save_sedflux_grids && push!(sel, :sedflux)
    p.save_sandfrac_grids && push!(sel, :sandfrac)
    if p.save_discharge_components
        push!(sel, :discharge_x, :discharge_y)
    end
    if p.save_velocity_components
        push!(sel, :velocity_x, :velocity_y)
    end
    return sel
end

"""Whether any NetCDF output is requested (grids always bring the metadata with them)."""
writes_netcdf(p::DeltaParameters) = p.save_metadata || !isempty(saved_grids(p))

"""
Auxiliary variables written once (`:scalar`, `:grid`) or at every save (`:series`), as
`(name, units, type, kind, long_name)`; `value(model, name)` supplies the value.
"""
const META_OUTPUTS = (
    (:L0, "cells", Int64, :scalar, "channel_entrance__length"),
    (:N0, "cells", Int64, :scalar, "channel_entrance__width"),
    (:CTR, "cells", Int64, :scalar, "channel_entrance__y_position"),
    (:dx, "meters", Float32, :scalar, "model_grid_cell_edge__length"),
    (:h0, "meters", Float32, :scalar, "channel_entrance__depth"),
    (:hb, "meters", Float32, :scalar, "basin_bottom_initial__depth"),
    (:cell_type, "type", Int64, :grid, "model_grid_cell__type"),
    (:start_subsidence, "seconds", Int64, :scalar, "basin_bottom_vertical_rate_of_change__start_time"),
    (:sigma, "meters per timestep", Float32, :grid, "basin_bottom__vertical_rate_of_change"),
    (:H_SL, "meters", Float32, :series, "basin_water_surface__elevation"),
    (:f_bedload, "fraction", Float32, :series, "channel_entrance_water_sediment_sand__volume_fraction"),
    (:C0_percent, "percent", Float32, :series, "channel_entrance__water_sediment__volume_percent"),
    (:u0, "meters per second", Float32, :series, "channel_entrance__speed"),
)

"""
    metadata_outputs(params)

The `META_OUTPUTS` entries applicable to a configuration (subsidence entries only when
subsidence is on).
"""
metadata_outputs(p::DeltaParameters) =
    filter(v -> p.toggle_subsidence || v[1] ∉ (:start_subsidence, :sigma), collect(META_OUTPUTS))

# `CTR` is stored 0-based in files (as pyDeltaRCM writes it); in memory it is 1-based.
function metadata_value(m, name::Symbol)
    name === :CTR && return m.CTR - 1
    name === :start_subsidence && return round(Int, m.start_subsidence)
    return getproperty(m, name)
end

"""Fields that can be saved as figures (`save_*_figs`): `name => model field`."""
const FIGURE_OUTPUTS = (eta = :eta, stage = :stage, depth = :depth, discharge = :qw, velocity = :uw,
                        sedflux = :qs, sandfrac = :sand_frac)

"""Names of the fields selected by the `save_*_figs` options."""
function saved_figures(p::DeltaParameters)
    flags = (eta = p.save_eta_figs, stage = p.save_stage_figs, depth = p.save_depth_figs,
             discharge = p.save_discharge_figs, velocity = p.save_velocity_figs,
             sedflux = p.save_sedflux_figs, sandfrac = p.save_sandfrac_figs)
    return Symbol[k for k in keys(flags) if flags[k]]
end

# ---- extension entry points (methods are added by the extensions) -------------------------------

"""
    plot_field(model, field = :eta; kwargs...) -> Figure

Heatmap of a model field (`:eta`, `:depth`, `:uw`, …) with the inlet at the top. Provided by
the Makie extension: load a backend (`using CairoMakie`) first.
"""
function plot_field end
"""Write the figure of `field` for output record `idx` into `out_dir`. Provided by the Makie extension."""
function save_figure end

"""Open (or, when resuming, reopen) the NetCDF output file. Provided by the NCDatasets extension."""
function open_netcdf! end
"""Append the current state to the NetCDF output. Provided by the NCDatasets extension."""
function write_netcdf! end
"""Close the NetCDF output file. Provided by the NCDatasets extension."""
function close_netcdf! end
"""
    save_checkpoint(model, path)

Write a restartable checkpoint. Provided by the JLD2 extension (`using JLD2`).
"""
function save_checkpoint end
"""
    load_checkpoint!(model, path)

Restore state written by [`save_checkpoint`](@ref). Provided by the JLD2 extension.
"""
function load_checkpoint! end

function _require(f, pkg, what)
    isempty(methods(f)) && throw(ErrorException("$what needs the $pkg package: run `using $pkg` first"))
    return nothing
end

# ---- driving the output from the model ----------------------------------------------------------

checkpoint_file(m) = joinpath(m.params.out_dir, "checkpoint.jld2")
output_file(m) = joinpath(m.params.out_dir, "pyDeltaRCM_output.nc")
checkpoint_interval(p::DeltaParameters) = something(p.checkpoint_dt, p.save_dt)

"""
Prepare output for a new or resumed model: create the output directory/file, and (for a
new model) record the initial state.
"""
function init_output!(m; resumed::Bool = false)
    p = m.params
    figs = saved_figures(p)
    if writes_netcdf(p) || p.save_checkpoint || !isempty(figs)
        mkpath(p.out_dir)
    end
    isempty(figs) || _require(save_figure, "Makie (e.g. CairoMakie)", "Figure output (`save_*_figs`)")
    if writes_netcdf(p)
        _require(open_netcdf!, "NCDatasets", "NetCDF output (`save_*_grids`, `save_metadata`)")
        m.output = open_netcdf!(m; append = resumed)
    end
    if !resumed
        m.since_data = Inf              # force a save at t = 0
        m.since_checkpoint = Inf
        output_data!(m)
        output_checkpoint!(m)
    end
    return m
end

function output_data!(m)
    m.since_data >= m.params.save_dt || return nothing
    if m.output !== nothing
        write_netcdf!(m.output, m, m.save_iter)
    end
    for name in saved_figures(m.params)
        save_figure(m, getproperty(FIGURE_OUTPUTS, name), m.save_iter)
    end
    m.save_iter += 1
    m.since_data = 0.0
    return nothing
end

function output_checkpoint!(m)
    m.params.save_checkpoint || return nothing
    m.since_checkpoint >= checkpoint_interval(m.params) || return nothing
    _require(save_checkpoint, "JLD2", "Checkpointing (`save_checkpoint`)")
    save_checkpoint(m, checkpoint_file(m))
    if checkpoint_interval(m.params) != m.params.save_dt && writes_netcdf(m.params)
        @warn "Grid save interval and checkpoint interval differ: resuming may duplicate entries in the NetCDF output"
    end
    m.since_checkpoint = 0.0
    return nothing
end

function close_output!(m)
    if m.output !== nothing
        close_netcdf!(m.output)
        m.output = nothing
    end
    return nothing
end

"""
    resume!(model, path = <out_dir>)

Restore `model` from the checkpoint in directory `path` (done automatically when the
`resume_checkpoint` parameter is set).
"""
function resume!(m, dir::AbstractString = m.params.out_dir)
    _require(load_checkpoint!, "JLD2", "Resuming from a checkpoint")
    file = joinpath(dir, "checkpoint.jld2")
    isfile(file) || throw(ArgumentError("No checkpoint found at $file"))
    load_checkpoint!(m, file)
    return m
end
