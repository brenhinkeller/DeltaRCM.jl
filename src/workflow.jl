# Running models for a duration, from configuration files, and in ensembles/parameter sweeps.

"""
    run!(model; timesteps, time, time_years, If = 1.0)

Run until a duration has elapsed (exactly one of the keywords): `timesteps` whole time steps,
`time` model seconds, or `time_years` real years. `time` and `time_years` are scaled by the
intermittency factor `If` (model time = real time × `If`); any remainder shorter than `dt` is
taken as a final shortened step, so the run ends exactly at the requested time.
"""
function run!(m::DeltaModel; timesteps = nothing, time = nothing, time_years = nothing, If::Real = 1.0)
    given = count(!isnothing, (timesteps, time, time_years))
    given == 1 || throw(ArgumentError("give exactly one of `timesteps`, `time`, `time_years`"))
    timesteps !== nothing && return run!(m, timesteps)
    duration = time !== nothing ? time * If : time_years * If * 86400 * 365.25
    duration >= 0 || throw(ArgumentError("run duration must be non-negative"))
    target = m.time + duration
    dt = m.c.dt
    nwhole = floor(Int, duration / dt * (1 + 1e-12))
    run!(m, nwhole)
    remainder = target - m.time
    remainder > 1e-9 * dt && step!(m; dt = remainder)
    return m
end

const DURATION_KEYS = (:timesteps, :time, :time_years)

"""
    run_model(config = nothing; kwargs...)

Configure, run and finalize one model, like `pyDeltaRCM`'s `Preprocessor`. `config` is a
YAML path or a `Dict`; keywords add to or override it. Besides the [`DeltaParameters`](@ref)
fields, it accepts the duration keys of [`run!`](@ref) (`timesteps`, `time`, `time_years`,
`If`), and `T` (the float type, default `Float64`). When the configuration has
`matrix`, `set` or `ensemble` entries, it is expanded to several jobs ([`run_jobs`](@ref)).
Returns the (finalized) model, or for several jobs the vector of results.
"""
function run_model(config = nothing; T::Type{<:AbstractFloat} = Float64, kwargs...)
    cfg = _config_dict(config; kwargs...)
    any(k -> haskey(cfg, k), (:matrix, :set, :ensemble)) && return run_jobs(identity, cfg; T)
    return _run_job(identity, cfg, T)
end

function _config_dict(config; kwargs...)
    cfg = config === nothing ? Dict{Symbol,Any}() :
          config isa AbstractString ? _load_yaml(config) :
          Dict{Symbol,Any}(Symbol(k) => v for (k, v) in config)
    for (k, v) in kwargs
        cfg[k] = v
    end
    return cfg
end

# one job: build, run for the configured duration, apply `f(model)`, finalize
function _run_job(f, cfg::Dict{Symbol,Any}, T)
    cfg = copy(cfg)
    duration = Dict{Symbol,Any}(k => pop!(cfg, k) for k in (DURATION_KEYS..., :If) if haskey(cfg, k))
    any(k -> haskey(duration, k), DURATION_KEYS) ||
        throw(ArgumentError("a run duration is required: one of `timesteps`, `time`, `time_years`"))
    foreach(k -> pop!(cfg, k, nothing), (:matrix, :set, :ensemble, :config, :dryrun, :parallel))
    m = DeltaModel{T}(DeltaParameters(cfg))
    try
        run!(m; (k => v for (k, v) in duration)...)
        return f(m)
    finally
        m.is_finalized || DeltaRCM.finalize!(m)
    end
end

# ---- multi-job expansion: matrix, set, ensemble -----------------------------------------------------

"""
    expand_jobs(config; rng = Random.default_rng()) -> Vector{Dict{Symbol,Any}}

Expand the `matrix` (Cartesian product), `set` (explicit list) and `ensemble` (replicates
with random `seed`s) entries of a configuration, as `pyDeltaRCM`'s preprocessor does, into
one configuration per job. With several jobs, `out_dir` is required and each job writes to
`out_dir/job_000`, `out_dir/job_001`, ….
"""
function expand_jobs(config::AbstractDict; rng::AbstractRNG = Random.default_rng())
    cfg = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in config)
    has_set, has_matrix, has_ens = haskey(cfg, :set), haskey(cfg, :matrix), haskey(cfg, :ensemble)
    (has_set && has_matrix) && throw(ArgumentError("Cannot use \"matrix\" and \"set\" expansion together"))
    if !(has_set || has_matrix || has_ens)
        return [cfg]
    end
    haskey(cfg, :out_dir) || throw(ArgumentError("You must specify `out_dir` to use any multi-job expansion"))
    nens = has_ens ? pop!(cfg, :ensemble) : 1
    nens isa Integer || throw(ArgumentError("Invalid ensemble type, must be an integer"))
    newseed() = Int(rand(rng, UInt32))
    combos = Dict{Symbol,Any}[]
    if has_set
        sets = pop!(cfg, :set)
        sets isa AbstractVector && all(s -> s isa AbstractDict, sets) ||
            throw(ArgumentError("`set` must be a list of mappings"))
        sets = [Dict{Symbol,Any}(Symbol(k) => v for (k, v) in s) for s in sets]
        all(s -> keys(s) == keys(first(sets)), sets) || throw(ArgumentError("All keys in all sets must be identical"))
        for s in sets, _ in 1:(nens > 1 ? nens : 1)
            push!(combos, nens > 1 ? merge(s, Dict(:seed => newseed())) : s)
        end
    else
        matrix = has_matrix ? Dict{Symbol,Any}(Symbol(k) => v for (k, v) in pop!(cfg, :matrix)) : Dict{Symbol,Any}()
        if nens > 1
            haskey(matrix, :seed) && throw(ArgumentError("Random seeds cannot be in the matrix when `ensemble` is given"))
            matrix[:seed] = [newseed() for _ in 1:nens]
        end
        haskey(matrix, :out_dir) && throw(ArgumentError("`out_dir` cannot be a matrix key"))
        for k in keys(matrix)
            haskey(cfg, k) && throw(ArgumentError("Key `$k` is in both the matrix and the fixed configuration"))
        end
        names = collect(keys(matrix))
        # first matrix key varies slowest, matching itertools.product in pyDeltaRCM
        for vals in Iterators.product(reverse([matrix[k] for k in names])...)
            push!(combos, Dict{Symbol,Any}(zip(reverse(names), vals)))
        end
    end
    root = cfg[:out_dir]
    return map(enumerate(combos)) do (i, combo)
        merge(cfg, combo, Dict(:out_dir => joinpath(root, "job_" * lpad(i - 1, 3, '0'))))
    end
end

"""
    run_jobs([f,] config; parallel = false, T = Float64, rng) -> Vector

Expand `config` with [`expand_jobs`](@ref) and run every job (each also finalized).
`f(model)` is applied to each model after its run, before finalization (e.g. to extract a
result); the returned vector holds `f`'s results (the models themselves by default).
With `parallel = true` jobs run as concurrent tasks (each job itself is single-threaded). A table of job parameters is written to
`jobs_parameters.txt` in the root `out_dir`.
"""
function run_jobs(f, config; parallel::Bool = false, T::Type{<:AbstractFloat} = Float64,
                  rng::AbstractRNG = Random.default_rng(), kwargs...)
    cfg = _config_dict(config; kwargs...)
    jobs = expand_jobs(cfg; rng)
    if any(k -> haskey(cfg, k), (:matrix, :set, :ensemble))
        root = cfg[:out_dir]
        if isdir(root) && get(cfg, :resume_checkpoint, false) !== true
            throw(ArgumentError("Job output directory ($root) already exists"))
        end
        mkpath(root)
        keys_varied = sort!(collect(filter(k -> any(j -> get(j, k, nothing) != get(first(jobs), k, nothing), jobs), keys(first(jobs)))))
        filter!(!=(:out_dir), keys_varied)
        open(joinpath(root, "jobs_parameters.txt"), "w") do io
            println(io, join(["job_id"; String.(keys_varied)], ", "))
            for j in jobs
                println(io, join([basename(j[:out_dir]); [string(get(j, k, "")) for k in keys_varied]], ","))
            end
        end
    end
    results = Vector{Any}(undef, length(jobs))
    if parallel
        @sync for (i, j) in enumerate(jobs)
            Threads.@spawn results[i] = _run_job(f, j, T)
        end
    else
        for (i, j) in enumerate(jobs)
            results[i] = _run_job(f, j, T)
        end
    end
    return results
end
run_jobs(config; kwargs...) = run_jobs(identity, config; kwargs...)
