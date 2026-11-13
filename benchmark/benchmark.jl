# Per-step timing of the default 100 × 200 domain (compare with `pydeltarcm_timing.py`).
#   julia --project=. benchmark/benchmark.jl
using DeltaRCM

m = DeltaModel(; seed = 0)
run!(m, 3)                                  # warm up
t = @elapsed run!(m, 20)
a = @allocated step!(m)
println(round(1000t / 20; digits = 1), " ms/step, ", a, " bytes allocated/step")
