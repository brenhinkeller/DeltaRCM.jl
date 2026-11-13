"""
    DeltaParameters

User-facing model configuration. Field names and defaults follow `pyDeltaRCM`'s
`default.yml` exactly, so that existing YAML configuration files work unchanged.

Construct from keywords, a YAML file, or a `Dict`:

```julia
DeltaParameters(; Np_water = 1000, seed = 0)
DeltaParameters("config.yaml"; seed = 1)       # keywords override the file
```
"""
Base.@kwdef struct DeltaParameters
    out_dir::String = "deltaRCM_Output"
    verbose::Int = 0
    seed::Union{Int,Nothing} = nothing
    Length::Float64 = 5000.0
    Width::Float64 = 10000.0
    dx::Float64 = 50.0
    L0_meters::Float64 = 150.0
    S0::Float64 = 0.00015
    itermax::Int = 3
    Np_water::Int = 2000
    u0::Float64 = 1.0
    N0_meters::Float64 = 250.0
    h0::Float64 = 5.0
    hb::Union{Float64,Nothing} = nothing
    H_SL::Float64 = 0.0
    SLR::Float64 = 0.0
    Np_sed::Int = 2000
    f_bedload::Float64 = 0.5
    active_layer_thickness::Union{Float64,Nothing} = nothing
    C0_percent::Float64 = 0.1
    Csmooth::Float64 = 0.9
    toggle_subsidence::Bool = false
    subsidence_rate::Float64 = 2e-9
    start_subsidence::Float64 = 216000.0
    save_eta_figs::Bool = false
    save_stage_figs::Bool = false
    save_depth_figs::Bool = false
    save_discharge_figs::Bool = false
    save_velocity_figs::Bool = false
    save_sedflux_figs::Bool = false
    save_sandfrac_figs::Bool = false
    save_figs_sequential::Bool = true
    save_metadata::Bool = false
    save_eta_grids::Bool = false
    save_stage_grids::Bool = false
    save_depth_grids::Bool = false
    save_discharge_grids::Bool = false
    save_velocity_grids::Bool = false
    save_sedflux_grids::Bool = false
    save_sandfrac_grids::Bool = false
    save_discharge_components::Bool = false
    save_velocity_components::Bool = false
    save_dt::Float64 = 86400.0
    checkpoint_dt::Union{Float64,Nothing} = nothing
    save_checkpoint::Bool = false
    resume_checkpoint::Union{Bool,String} = false
    omega_sfc::Float64 = 0.1
    omega_flow::Float64 = 0.9
    Nsmooth::Int = 10
    theta_water::Float64 = 1.0
    coeff_theta_sand::Float64 = 2.0
    coeff_theta_mud::Float64 = 1.0
    beta::Int = 3
    sed_lag::Float64 = 1.0
    coeff_U_dep_mud::Float64 = 0.3
    coeff_U_ero_mud::Float64 = 1.5
    coeff_U_ero_sand::Float64 = 1.05
    alpha::Float64 = 0.1
    stepmax::Union{Int,Nothing} = nothing
    force_deposit::Bool = false
    sand_frac_bc::Float64 = 0.0
    clobber_netcdf::Bool = false
    legacy_netcdf::Bool = false
end

# Keys that pyDeltaRCM's preprocessor consumes; ignored here without warning.
const PREPROCESSOR_KEYS = (:timesteps, :time, :time_years, :config, :dryrun, :parallel)

# ---- conversion of user values to field types -------------------------------------------

_convert_param(::Type{T}, name, v) where {T<:AbstractString} = String(v)
function _convert_param(::Type{Bool}, name, v)
    v isa Bool || throw(ArgumentError("Parameter `$name` must be a Bool, got $(repr(v))"))
    return v
end
function _convert_param(::Type{Int}, name, v)
    (v isa Real && !(v isa Bool) && isinteger(v)) ||
        throw(ArgumentError("Parameter `$name` must be an integer, got $(repr(v))"))
    return Int(v)
end
function _convert_param(::Type{Float64}, name, v)
    (v isa Real && !(v isa Bool)) ||
        throw(ArgumentError("Parameter `$name` must be a number, got $(repr(v))"))
    return Float64(v)
end
function _convert_param(::Type{Union{A,Nothing}}, name, v) where {A}
    (v === nothing || v === missing || v == "None" || v == "null") && return nothing
    return _convert_param(A, name, v)
end
function _convert_param(::Type{Union{Bool,String}}, name, v)
    v isa Bool && return v
    v isa AbstractString && return String(v)
    throw(ArgumentError("Parameter `$name` must be a Bool or String, got $(repr(v))"))
end

function _validate(p::DeltaParameters)
    p.dx > 0 || throw(ArgumentError("`dx` must be positive"))
    p.Length > 0 && p.Width > 0 || throw(ArgumentError("`Length` and `Width` must be positive"))
    p.itermax >= 1 || throw(ArgumentError("`itermax` must be at least 1"))
    p.Np_water >= 1 || throw(ArgumentError("`Np_water` must be at least 1"))
    p.Np_sed >= 1 || throw(ArgumentError("`Np_sed` must be at least 1"))
    p.u0 > 0 && p.h0 > 0 || throw(ArgumentError("`u0` and `h0` must be positive"))
    0 <= p.f_bedload <= 1 || throw(ArgumentError("`f_bedload` must be in [0, 1]"))
    0 <= p.omega_sfc <= 1 || throw(ArgumentError("`omega_sfc` must be in [0, 1]"))
    0 <= p.omega_flow <= 1 || throw(ArgumentError("`omega_flow` must be in [0, 1]"))
    0 <= p.Csmooth <= 1 || throw(ArgumentError("`Csmooth` must be in [0, 1]"))
    p.Nsmooth >= 0 || throw(ArgumentError("`Nsmooth` must be non-negative"))
    p.C0_percent > 0 || throw(ArgumentError("`C0_percent` must be positive"))
    p.S0 >= 0 || throw(ArgumentError("`S0` must be non-negative"))
    p.beta >= 1 || throw(ArgumentError("`beta` must be at least 1"))
    if p.active_layer_thickness !== nothing && p.active_layer_thickness < 0
        throw(ArgumentError("`active_layer_thickness` must be non-negative"))
    end
    if p.stepmax !== nothing && p.stepmax < 1
        throw(ArgumentError("`stepmax` must be at least 1"))
    end
    return p
end

"""
    DeltaParameters(config::AbstractDict; kwargs...)
    DeltaParameters(path::AbstractString; kwargs...)

Build parameters from a dictionary or YAML file; `kwargs` take precedence over the file
(with a warning, as in pyDeltaRCM). Unknown keys produce a warning.
"""
function DeltaParameters(config::AbstractDict; kwargs...)
    cfg = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in config)
    for (k, v) in kwargs
        haskey(cfg, k) && @warn "Keyword `$k` also specified in the configuration file; using the keyword value."
        cfg[k] = v
    end
    names = fieldnames(DeltaParameters)
    types = fieldtypes(DeltaParameters)
    vals = Dict{Symbol,Any}()
    unused = Symbol[]
    for (k, v) in cfg
        i = findfirst(==(k), names)
        if i === nothing
            k in PREPROCESSOR_KEYS || push!(unused, k)
        else
            vals[k] = _convert_param(types[i], k, v)
        end
    end
    isempty(unused) || @warn "Unused configuration keys: $(sort!(unused))"
    return _validate(DeltaParameters(; vals...))
end
DeltaParameters(path::AbstractString; kwargs...) =
    DeltaParameters(_load_yaml(path); kwargs...)

function _load_yaml(path::AbstractString)
    isfile(path) || throw(ArgumentError("Configuration file not found: $path"))
    d = YAML.load_file(path; dicttype = Dict{Symbol,Any})
    d === nothing && return Dict{Symbol,Any}()
    d isa AbstractDict || throw(ArgumentError("Configuration file $path is not a YAML mapping"))
    return d
end

Base.NamedTuple(p::DeltaParameters) =
    NamedTuple{fieldnames(DeltaParameters)}(ntuple(i -> getfield(p, i), fieldcount(DeltaParameters)))

"""
    with(p::DeltaParameters; kwargs...)

Return a copy of `p` with the given fields replaced.
"""
with(p::DeltaParameters; kwargs...) = _validate(DeltaParameters(; merge(NamedTuple(p), NamedTuple(kwargs))...))

function Base.show(io::IO, ::MIME"text/plain", p::DeltaParameters)
    d = DeltaParameters()
    println(io, "DeltaParameters (non-default values shown)")
    shown = false
    for n in fieldnames(DeltaParameters)
        v = getfield(p, n)
        if v != getfield(d, n)
            println(io, "  ", n, " = ", repr(v))
            shown = true
        end
    end
    shown || print(io, "  (all defaults)")
end
