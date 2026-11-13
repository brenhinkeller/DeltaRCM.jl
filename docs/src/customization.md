# Customization

pyDeltaRCM is customized by subclassing and overriding `hook_*` methods. Here callbacks are passed as
`hooks`, and are called with the model at these points: `:init`, `:before_step`, `:before_water`, `:after_water`, `:before_sediment`, `:after_sediment`, `:after_step`.

```julia
# make the left half of the basin less attractive to water, as a "tilted basin" experiment
model = DeltaModel(; Length = 2000, Width = 4000, dx = 50,
                   hooks = (; init = m -> (m.mod_water_weight[:, 1:m.CTR] .= 0.2),
                              after_step = m -> m.iter % 10 == 0 && @info "step" m.iter))
```

The arrays `mod_water_weight`, `mod_sed_weight`, `mod_erosion` and `mod_stable_weight` (all ones by default)
scale the corresponding weights cell by cell. For anything deeper, the building blocks
(`route_water!`, `route_sediment!`, `finalize_timestep!`, …) are plain functions on the model, and the
sediment and water routing dispatch on parcel kind, so new behaviors can be added as methods.
