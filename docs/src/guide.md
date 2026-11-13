# Guide

## Running a model

```julia
using DeltaRCM

model = DeltaModel(; Length = 2000, Width = 4000, dx = 50, seed = 0)
step!(model)               # one time step
run!(model, 50)            # 50 time steps
run!(model; time_years = 1, If = 0.1)   # until a model duration has elapsed
model.eta                  # bed elevation, an (L, W) array; first index is downstream
```

All fields (`eta`, `stage`, `depth`, `qx`, `qy`, `qw`, `ux`, `uy`, `uw`, `qs`, `sand_frac`, `cell_type`, …) are
`(L, W)` views of the model state, with the inlet at row 1. Parameters and derived constants are also
properties: `model.dt`, `model.Qs0`, `model.u0`, …

## Configuration

Parameters are those of pyDeltaRCM, with the same names and defaults ([`DeltaParameters`](@ref)). Pass them as
keywords, in a `Dict`, or in a YAML file:

```julia
DeltaModel("config.yaml"; seed = 3)        # keywords override the file (with a warning)
run_model("config.yaml")                   # build, run for `timesteps`/`time`/`time_years`, finalize
```

As in pyDeltaRCM, `Length`, `Width`, `L0_meters` and `N0_meters` are rounded to whole cells (with a warning).

## Precision

`DeltaModel{Float32}(...)` uses single precision. A model runs on a single thread; every parcel has its own
random stream derived from `(seed, time step, parcel)`, so results depend only on `seed`. Independent jobs
can run concurrently with `run_jobs(...; parallel = true)` (start Julia with `-t auto`).

## Ensembles and parameter sweeps

`matrix` (Cartesian product), `set` (explicit list) and `ensemble` (random seeds) work as in pyDeltaRCM, in a YAML
file or a `Dict`:

```julia
results = run_jobs(Dict(:out_dir => "sweep", :timesteps => 100, :Length => 2000, :Width => 4000,
                        :matrix => Dict(:f_bedload => [0.3, 0.5, 0.7]), :ensemble => 4)) do model
    sum(model.eta)            # applied to each finished model, before finalization
end
```

Each job writes to `out_dir/job_000`, `job_001`, …, and `jobs_parameters.txt` tabulates the varied parameters.
See [`expand_jobs`](@ref), [`run_jobs`](@ref).
