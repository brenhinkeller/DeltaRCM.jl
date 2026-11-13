module DeltaRCMNCDatasetsExt

using DeltaRCM
using NCDatasets
import DeltaRCM: open_netcdf!, write_netcdf!, close_netcdf!, NetCDFOutput

# File layout: pyDeltaRCM's "sandsuet 1.0.0" schema. Dimensions are (seconds, x, y) as seen by
# NetCDF/Python/xarray, with x downstream (`L`) and y cross-stream (`W`). NCDatasets lists
# dimensions in Julia (column-major) order, i.e. reversed, so a stored grid is `(W, L, time)`.

const SANDSUET_VERSION = "1.0.0"

timename(p) = p.legacy_netcdf ? "time" : "seconds"
groupname(p) = p.legacy_netcdf ? "meta" : "auxdata"
attribs(units, long) = ["units" => units, "long_name" => long]

function DeltaRCM.open_netcdf!(m; append::Bool = false)
    p = m.params
    path = DeltaRCM.output_file(m)
    if append
        if isfile(path)
            return NetCDFOutput(path, NCDataset(path, "a"))
        end
        @warn "NetCDF output file not found, but was expected. Creating a new output file."
        m.save_iter = 0
    end
    if isfile(path)
        p.clobber_netcdf || throw(ErrorException("Existing NetCDF output file at $path (set `clobber_netcdf` to replace it)"))
        @warn "Replacing existing netCDF file"
        rm(path)
    end
    ds = NCDataset(path, "c"; format = :netcdf4)
    ds.attrib["description"] = "Output from DeltaRCM.jl"
    ds.attrib["source"] = "DeltaRCM.jl v$(pkgversion(DeltaRCM))"
    ds.attrib["history"] = "Created " * Libc.strftime(time())
    p.legacy_netcdf || (ds.attrib["sandsuet_version"] = SANDSUET_VERSION)

    L, W = size(m)
    tname = timename(p)
    defDim(ds, "x", L)
    defDim(ds, "y", W)
    defDim(ds, tname, Inf)
    defVar(ds, tname, Float32, (tname,); attrib = ["units" => "second"])
    x = defVar(ds, "x", Float32, ("x",); attrib = ["units" => "meter"])
    y = defVar(ds, "y", Float32, ("y",); attrib = ["units" => "meter"])
    x[:] = Float32.((0:L-1) .* m.dx)
    y[:] = Float32.((0:W-1) .* m.dx)

    for name in DeltaRCM.saved_grids(p)
        _, units, long = getproperty(DeltaRCM.GRID_OUTPUTS, name)
        defVar(ds, String(name), Float32, ("y", "x", tname); attrib = attribs(units, long))
    end

    grp = defGroup(ds, groupname(p))
    for (name, units, T, kind, long) in DeltaRCM.metadata_outputs(p)
        dims = kind === :scalar ? () : kind === :grid ? ("y", "x") : (tname,)
        v = defVar(grp, String(name), T, dims; attrib = attribs(units, long))
        value = DeltaRCM.metadata_value(m, name)
        if kind === :scalar
            v[] = T(value)
        elseif kind === :grid
            v[:, :] = T.(permutedims(value))
        end
    end
    return NetCDFOutput(path, ds)
end

function DeltaRCM.write_netcdf!(out::NetCDFOutput, m, idx::Integer)
    ds = out.dataset
    p = m.params
    k = idx + 1
    ds[timename(p)][k] = Float32(m.time)
    for name in DeltaRCM.saved_grids(p)
        field = first(getproperty(DeltaRCM.GRID_OUTPUTS, name))
        ds[String(name)][:, :, k] = Float32.(permutedims(getproperty(m, field)))
    end
    grp = ds.group[groupname(p)]
    for (name, _, T, kind, _) in DeltaRCM.metadata_outputs(p)
        kind === :series && (grp[String(name)][k] = T(DeltaRCM.metadata_value(m, name)))
    end
    NCDatasets.sync(ds)
    return nothing
end

DeltaRCM.close_netcdf!(out::NetCDFOutput) = (close(out.dataset); nothing)

end
