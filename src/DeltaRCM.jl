"""
    DeltaRCM

A reduced-complexity river-delta model (Liang et al. 2015) in pure Julia, with the user
interface and configuration of `pyDeltaRCM`.

```julia
using DeltaRCM
model = DeltaModel(; Length = 2000, Width = 4000, dx = 50, seed = 0)
run!(model, 50)
model.eta
```
"""
module DeltaRCM

using PrecompileTools
using Random
using YAML

export DeltaParameters, DeltaModel, step!, run!, finalize!, resume!, run_model, run_jobs, expand_jobs, plot_field

include("parameters.jl")
include("stencil.jl")
include("domain.jl")
include("water.jl")
include("sediment.jl")
include("io.jl")
include("model.jl")
include("workflow.jl")

# Compile the whole time step ahead of time (tiny domain), for both float types.
@setup_workload begin
    tiny = (; Length = 1000, Width = 2000, dx = 50, Np_water = 20, Np_sed = 20, seed = 0)
    @compile_workload begin
        for T in (Float64, Float32)
            m = DeltaModel{T}(; tiny...)
            run!(m, 2)
            step!(m; dt = m.c.dt / 2)
            finalize!(m)
        end
        show(devnull, DeltaModel(; tiny...))
        show(devnull, DeltaParameters())
    end
end

end
