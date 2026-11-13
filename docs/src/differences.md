# Differences from pyDeltaRCM

Same model, same parameters and file formats; different numerics where the existing codes carried
quirks that are not part of the model as published. Over ensembles of runs, statistics of the delta (shoreline
position, land area, flow speed and depth, surface sand fraction) agree with pyDeltaRCM; individual runs differ
because the random streams differ.

Deliberate differences:

- Every parcel has its own random stream, derived from `(seed, time step, parcel)`.
- The domain has a wall halo rather than edge padding; parcels reaching an `EDGE` cell are exported without
  depositing or eroding there, and the flux through edge cells is not double counted.
- The flooding correction is a two-phase calculation of wetness from the current stage and bed.
- Topographic diffusion excludes wall cells and the inlet row, so it conserves bed mass exactly.
- The water surface is always under-relaxed (as in the MATLAB model), including the first step.
- Velocity is zero wherever the depth does not exceed `dry_depth`.
- Mass conservation is exact: bed-volume change = sediment in − exported − lost + inlet boundary term.
- Loops in water parcel paths are detected in O(1) with visit stamps.

Kept from the published model: the `omega_flow` schedule, the loop-escape relocation, the `gamma`
stability check, and the weighting formulae.
