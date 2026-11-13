# DeltaRCM.jl

An attempt to accelerate reduced complexity delta models in Julia.

A fast, pure-Julia reduced-complexity river-delta model: the [DeltaRCM](https://github.com/csdms-contrib/DeltaRCM) of
Liang et al. ([2015a](https://doi.org/10.5194/esurf-3-67-2015), [2015b](https://doi.org/10.5194/esurf-3-87-2015)), 
with the configuration and output conventions of [pyDeltaRCM](https://github.com/DeltaRCM/pyDeltaRCM)
([Moodie et al. (2021)](https://doi.org/10.21105/joss.03398)).

The initial implementation is heavily LLM-assisted (Claude code, Sonnet 5.5),
and informed by both the original [DeltaRCM](https://github.com/csdms-contrib/DeltaRCM) (Matlab, GPL) 
and [pyDeltaRCM](https://github.com/DeltaRCM/pyDeltaRCM) (Python, MIT license) open-source implementations.
We retain the GPL of the original DeltaRCM implementation as a result.

```julia
using DeltaRCM

model = DeltaModel(; Length = 5000, Width = 10000, dx = 50, seed = 42)   # pyDeltaRCM parameter names
run!(model; timesteps = 100)
model.eta                                                                 # bed elevation, (L, W)
```

or from a pyDeltaRCM YAML file:

```julia
run_model("config.yaml")                   # `timesteps`/`time`/`time_years`, `matrix`, `set`, `ensemble`
```

see [examples] for more examples, including plotting.

## Features

- **Same interface as pyDeltaRCM.** The same parameters and YAML keys (`Np_water`, `f_bedload`, `save_eta_grids`, …), the same `matrix`/`set`/`ensemble` job expansion, `save_dt`, checkpoints and resuming.
- **Fast.** About 33 ms per time step on one thread for the default domain (100 x 200 grid, Apple M1 Max), with zero allocations per step. Ensembles and sweeps can run concurrently (`run_jobs(...; parallel = true)`).
- **Reproducible.** Every parcel draws from its own random stream derived from `(seed, time step, parcel)`, so results depend only on `seed`.
- **Pure Julia**, `Float64` or `Float32` (`DeltaModel{Float32}(...)`).
- **NetCDF output** in pyDeltaRCM's layout (`using NCDatasets`), readable by existing xarray/netCDF4 workflows; **checkpoints** with JLD2; **plots** and `save_*_figs` with Makie.
- **Customizable** through callbacks (`hooks = (; before_water = m -> ...)`) and the `mod_*` weighting arrays.

## Differences from pyDeltaRCM

Results agree statistically (shoreline position, land area, flow speed and depth, surface sand fraction over ensembles), but not trajectory-for-trajectory: the random numbers differ, and a handful of potentially unintentional quirks of the original implementations were deliberately not reproduced (see the documentation for the full list).


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://brenhinkeller.github.io/DeltaRCM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://brenhinkeller.github.io/DeltaRCM.jl/dev/)
[![Build Status](https://github.com/brenhinkeller/DeltaRCM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brenhinkeller/DeltaRCM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/brenhinkeller/DeltaRCM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/brenhinkeller/DeltaRCM.jl)
