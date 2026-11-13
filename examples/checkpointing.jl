# Checkpointing: interrupting a run and resuming it exactly.
#
#   julia --project=examples examples/checkpointing.jl
#
# (first run: `julia --project=examples -e 'using Pkg; Pkg.instantiate()'`)
#
# Writes to `examples/output/checkpoint_demo/`: `checkpoint.jld2` and `pyDeltaRCM_output.nc`.
#
# With `save_checkpoint = true`, the model rewrites `<out_dir>/checkpoint.jld2` every
# `checkpoint_dt` model seconds (default: every `save_dt`). Setting `resume_checkpoint = true`
# restores the state from that file. Because every random number in the model derives from
# `(seed, time step, parcel)`, no random-number state needs saving, and a resumed run is
# bit-for-bit identical to one that was never interrupted. Checkpoints need JLD2, and the
# NetCDF output used here needs NCDatasets: loading each package enables the feature.

using DeltaRCM
using JLD2
using NCDatasets

outdir = joinpath(@__DIR__, "output", "checkpoint_demo")
rm(outdir; force = true, recursive = true)

total_steps = 120

# Everything that defines the model. The seed must be given for resumed and reference runs to agree.
setup = (; Length = 3000, Width = 6000, dx = 50, Np_water = 1000, Np_sed = 1000, seed = 42)
# Grids and a checkpoint every `save_dt` = 1 model day, which is every 4th step (dt = 25 000 s).
io = (; out_dir = outdir, save_dt = 86400, save_eta_grids = true, save_checkpoint = true)

## --- 1. a reference run that is never interrupted -------------------------------------------
reference = DeltaModel(; setup...)             # no output needed
run!(reference; timesteps = total_steps)

## --- 2. a run that is interrupted -------------------------------------------------------------
first_run = DeltaModel(; setup..., io...)
run!(first_run; timesteps = 50)
println("interrupted at step ", first_run.iter, " (time ", round(first_run.time / 86400; digits = 2), " days)")
finalize!(first_run)                            # as if the job had stopped here

## --- 3. resume --------------------------------------------------------------------------------------
# The last checkpoint was written at step 48, so that is where we resume; steps 49 and 50 are
# simply computed again. `resume_checkpoint` may also be the path to another run's output directory.
resumed = DeltaModel(; setup..., io..., resume_checkpoint = true)
println("resumed at step ", resumed.iter, " (time ", round(resumed.time / 86400; digits = 2), " days)")
run!(resumed; timesteps = total_steps - resumed.iter)
finalize!(resumed)

## --- 4. check -------------------------------------------------------------------------------------------
println("\nafter ", total_steps, " steps, resumed run vs uninterrupted run:")
for name in (:eta, :stage, :depth, :qx, :qy, :uw, :sand_frac)
    same = getproperty(resumed, name) == getproperty(reference, name)
    println("  ", rpad(name, 10), same ? "identical" : "DIFFERENT")
end

# The NetCDF file continues seamlessly too: the resumed run reopens it and carries on from the
# record at the checkpoint, overwriting anything the interrupted run wrote after it.
NCDataset(joinpath(outdir, "pyDeltaRCM_output.nc")) do ds
    t = ds["seconds"][:] ./ 86400
    println("\nNetCDF records: ", length(t), " (every ", round(t[2] - t[1]; digits = 2), " days up to day ",
            round(t[end]; digits = 2), "), strictly increasing: ", issorted(t; lt = <=))
    println("final record matches the model: ", maximum(abs.(permutedims(ds["eta"][:, :, end]) .- reference.eta)) < 1e-4)
end

## --- End of File