# Output and checkpoints

## NetCDF

With any `save_*_grids` option or `save_metadata`, loading `NCDatasets` enables output to
`<out_dir>/pyDeltaRCM_output.nc` every `save_dt` model seconds (and at `t = 0`), in pyDeltaRCM's
"sandsuet 1.0.0" layout:

| variable | contents |
|---|---|
| `eta`, `stage`, `depth` | bed elevation, water surface, water depth (m) |
| `discharge`, `velocity`, `sedflux`, `sandfrac` | `qw`, `uw`, `qs`, `sand_frac` |
| `discharge_x/y`, `velocity_x/y` | components (with `save_discharge_components`, `save_velocity_components`) |
| `seconds`, `x`, `y` | coordinates (`time` instead of `seconds` with `legacy_netcdf`) |
| `auxdata/…` (`meta/…` if `legacy_netcdf`) | `L0`, `N0`, `CTR`, `dx`, `h0`, `hb`, `cell_type`, `H_SL`, … |

Grids have NetCDF dimensions `(seconds, x, y)` with `x` downstream, so `xarray.open_dataset(...).eta` works
as for pyDeltaRCM output. (`CTR` is written 0-based, as in pyDeltaRCM.)

## Checkpoints and resuming

With `save_checkpoint = true` (and `using JLD2`), `<out_dir>/checkpoint.jld2` is rewritten every `checkpoint_dt`
(default `save_dt`). Set `resume_checkpoint = true` (or the path of a directory containing a checkpoint) to
continue: the resumed run is **bit-identical** to an uninterrupted one, since random numbers depend only on
`seed` and time step. An existing NetCDF file is extended in place.

## Figures

`save_eta_figs` etc. write PNGs to `out_dir` (with `using CairoMakie`); [`plot_field`](@ref) makes a figure of
any field on demand.
