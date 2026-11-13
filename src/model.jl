"""
    DeltaModel{T}(config; hooks, kwargs...)

A delta model with state at precision `T` (default `Float64`). `config` may be a
[`DeltaParameters`](@ref), a path to a pyDeltaRCM-style YAML file, a `Dict`, or omitted
(then everything comes from keywords):

```julia
model = DeltaModel(; Length = 2000, Width = 4000, dx = 50, seed = 42)
model = DeltaModel("config.yaml"; Np_water = 500)
step!(model)            # one time step
run!(model, 100)        # 100 time steps
model.eta               # bed elevation, an (L, W) view (first index downstream)
```

Results are reproducible for a given `seed`: every parcel draws from its own random stream.
"""
mutable struct DeltaModel{T<:AbstractFloat}
    params::DeltaParameters
    c::Constants{T}
    f::Fields{T}
    inlet::Vector{Int}                 # padded column indices of the inlet cells (row 1)
    scratch::WaterScratch{T}
    flood::FloodBuffer{T}
    tally::SedimentTally{T}
    seed::UInt64
    time::Float64
    iter::Int
    H_SL::T
    Vp_inletbc::T                      # bed volume added at the inlet by its boundary condition last step
    is_finalized::Bool
    save_iter::Int                     # index of the next record in the output file
    since_data::Float64                # model seconds since grids were last saved
    since_checkpoint::Float64
    output::Union{Nothing,NetCDFOutput}
    dt0::T                             # the nominal time step (a final partial step may be shorter)
    hooks::Any                         # NamedTuple of callbacks `f(model)`, see `HOOK_POINTS`
end

function DeltaModel{T}(params::DeltaParameters; hooks = (;)) where {T<:AbstractFloat}
    _check_hooks(hooks)
    p = resolve_parameters(params)
    c = Constants{T}(p)
    f = Fields{T}(c.Lp, c.Wp)
    inlet = init_domain!(f, c)
    scratch = WaterScratch(f, c)
    seed = p.seed === nothing ? rand(Random.RandomDevice(), UInt64) : UInt64(p.seed)
    if p.toggle_subsidence
        for j in 1:c.W, i in (c.L0+1):c.L
            f.sigma[i+1, j+1] = c.subsidence_rate * c.dt
        end
    end
    m = DeltaModel{T}(p, c, f, inlet, scratch, FloodBuffer{T}(), SedimentTally{T}(), seed, 0.0, 0,
                      T(p.H_SL), zero(T), false, 0, Inf, Inf, nothing, c.dt, hooks)
    resumed = p.resume_checkpoint !== false
    if resumed
        resume!(m, p.resume_checkpoint === true ? p.out_dir : String(p.resume_checkpoint))
    end
    call_hook(m, :init)
    init_output!(m; resumed)
    return m
end

DeltaModel{T}(config::Union{AbstractString,AbstractDict}; hooks = (;), kwargs...) where {T} =
    DeltaModel{T}(DeltaParameters(config; kwargs...); hooks)
DeltaModel{T}(; hooks = (;), kwargs...) where {T} = DeltaModel{T}(DeltaParameters(; kwargs...); hooks)
DeltaModel(args...; kwargs...) = DeltaModel{Float64}(args...; kwargs...)

# ---- property access: gridded fields as interior views, parameters, constants -----------------

const CONSTANT_NAMES = fieldnames(Constants)
const PARAMETER_NAMES = fieldnames(DeltaParameters)

@inline interior(A::AbstractMatrix) = @inbounds view(A, 2:size(A, 1)-1, 2:size(A, 2)-1)

const MODEL_FIELDS = fieldnames(DeltaModel)

# Real struct fields (`seed`, `H_SL`, `time`, ...) take precedence over same-named parameters.
function Base.getproperty(m::DeltaModel, s::Symbol)
    if s in MODEL_FIELDS
        return getfield(m, s)
    elseif s in FIELD_NAMES
        return interior(getfield(getfield(m, :f), s))
    elseif s in CONSTANT_NAMES
        return getfield(getfield(m, :c), s)
    elseif s in PARAMETER_NAMES
        return getfield(getfield(m, :params), s)
    elseif s === :time_iter
        return getfield(m, :iter)
    else
        throw(ArgumentError("DeltaModel has no property `$s`"))
    end
end

Base.propertynames(::DeltaModel) =
    (fieldnames(DeltaModel)..., FIELD_NAMES..., CONSTANT_NAMES..., PARAMETER_NAMES..., :time_iter)

Base.size(m::DeltaModel) = (m.c.L, m.c.W)

function Base.show(io::IO, m::DeltaModel{T}) where {T}
    print(io, "DeltaModel{$T}(", m.c.L, "×", m.c.W, " cells, dx=", m.c.dx, " m, t=", round(m.time; digits = 1),
          " s, step ", m.iter, ")")
end

# ---- customization hooks -----------------------------------------------------------------------

"""
Points at which a callback `f(model)` passed as `hooks = (; name = f, ...)` is called:
`:init` (after construction, before the initial output), `:before_step`, `:before_water`,
`:after_water`, `:before_sediment`, `:after_sediment` and `:after_step` (counters advanced,
before output is written). These correspond to pyDeltaRCM's `hook_*` methods. Hooks may read
and modify the model's fields, notably the weighting arrays `mod_water_weight`,
`mod_sed_weight`, `mod_erosion` and `mod_stable_weight` (all `1` by default).
"""
const HOOK_POINTS = (:init, :before_step, :before_water, :after_water, :before_sediment, :after_sediment,
                     :after_step)

function _check_hooks(hooks)
    hooks isa Union{NamedTuple,AbstractDict} || throw(ArgumentError("`hooks` must be a NamedTuple of functions"))
    for k in keys(hooks)
        Symbol(k) in HOOK_POINTS || throw(ArgumentError("Unknown hook `$k`; available: $(HOOK_POINTS)"))
    end
    return nothing
end

@inline function call_hook(m::DeltaModel, name::Symbol)
    h = getfield(m, :hooks)
    f = h isa NamedTuple ? get(h, name, nothing) : get(h, name, nothing)
    f === nothing || f(m)
    return nothing
end

# ---- the time step ------------------------------------------------------------------------

"""
    route_water!(model)

`itermax` iterations of: parcel routing, free-surface construction, flow update.
"""
function route_water!(m::DeltaModel)
    c, f = getfield(m, :c), getfield(m, :f)
    for iteration in 1:c.itermax
        run_water_parcels!(m.scratch, f, c, m.inlet, m.seed, m.iter, iteration, m.H_SL)
        compute_free_surface!(f, c, m.H_SL, m.flood)
        finalize_water_iteration!(f, c, m.inlet, iteration, m.iter, m.H_SL)
    end
    return m
end

"""
    route_sediment!(model)

Sand parcels, topographic diffusion, then mud parcels.
"""
function route_sediment!(m::DeltaModel)
    c, f = getfield(m, :c), getfield(m, :f)
    fill!(f.qs, 0); fill!(f.Vp_dep_sand, 0); fill!(f.Vp_dep_mud, 0)
    reset!(m.tally)
    route_sediment_parcels!(Sand(), f, c, m.inlet, m.seed, m.iter, m.tally)
    topo_diffusion!(f, c)
    route_sediment_parcels!(Mud(), f, c, m.inlet, m.seed, m.iter, m.tally)
    return m
end

function apply_subsidence!(m::DeltaModel)
    c, f = getfield(m, :c), getfield(m, :f)
    if m.toggle_subsidence && m.time >= c.start_subsidence
        f.eta .-= f.sigma .* (c.dt / m.dt0)
    end
    return m
end

"""Surface sand fraction as a continuously updated record of the bed (output only)."""
function compute_sand_frac!(f::Fields{T}, c::Constants{T}) where {T}
    a = c.active_layer_thickness
    @inbounds for j in 2:c.Wp-1, i in 2:c.Lp-1
        eta, eta0 = f.eta[i, j], f.eta0[i, j]
        if eta < eta0
            f.sand_frac[i, j] = eta < eta0 - a ? c.sand_frac_bc : f.active_layer[i, j]
        end
        dep = f.Vp_dep_sand[i, j] + f.Vp_dep_mud[i, j]
        if dep > 0
            mix = f.Vp_dep_sand[i, j] / dep
            f.sand_frac[i, j] = mix
            f.active_layer[i, j] = mix
        end
    end
    return nothing
end

function finalize_timestep!(m::DeltaModel{T}) where {T}
    c, f = getfield(m, :c), getfield(m, :f)
    flooding_correction!(f, c, m.flood)
    @inbounds for j in 2:c.Wp-1, i in 2:c.Lp-1
        s = max(f.stage[i, j], m.H_SL)
        f.stage[i, j] = s
        f.depth[i, j] = max(s - f.eta[i, j], zero(T))
    end
    # bed-elevation boundary condition at the inlet: bed sits h0 below the water surface
    Vbc = zero(T)
    @inbounds for jp in m.inlet
        Δ = (f.stage[2, jp] - c.h0) - f.eta[2, jp]
        Vbc += Δ * c.dx^2
        f.eta[2, jp] += Δ
        f.depth[2, jp] = c.h0
    end
    m.Vp_inletbc = Vbc
    compute_sand_frac!(f, c)
    m.H_SL += c.SLR * T(c.dt)
    return m
end

"""
    step!(model; dt = nothing)

Advance the model one time step: water routing, sediment routing, subsidence, and
boundary-condition bookkeeping. Giving `dt` (seconds) takes a step of that duration instead
of the nominal one (sediment volumes, diffusion and subsidence scale with it).
"""
function step!(m::DeltaModel; dt = nothing)
    (dt === nothing || dt == m.c.dt) && return _step!(m)
    c0 = getfield(m, :c)
    m.c = with_dt(c0, dt)
    try
        return _step!(m)
    finally
        m.c = c0
    end
end

function _step!(m::DeltaModel)
    m.is_finalized && throw(ErrorException("Cannot step a finalized model"))
    f = getfield(m, :f)
    call_hook(m, :before_step)
    copyto!(f.eta0, f.eta)
    call_hook(m, :before_water)
    route_water!(m)
    call_hook(m, :after_water)
    call_hook(m, :before_sediment)
    route_sediment!(m)
    call_hook(m, :after_sediment)
    apply_subsidence!(m)
    finalize_timestep!(m)
    m.time += m.c.dt
    m.iter += 1
    m.since_data += m.c.dt
    m.since_checkpoint += m.c.dt
    call_hook(m, :after_step)
    output_data!(m)
    output_checkpoint!(m)
    return m
end

"""
    run!(model, nsteps)

Take `nsteps` time steps.
"""
function run!(m::DeltaModel, nsteps::Integer)
    for _ in 1:nsteps
        step!(m)
    end
    return m
end

"""
    finalize!(model)

Close the output file. A finalized model cannot be stepped.
"""
function finalize!(m::DeltaModel)
    m.is_finalized && throw(ErrorException("Model already finalized"))
    close_output!(m)
    m.is_finalized = true
    return m
end
