module DeltaRCMJLD2Ext

using DeltaRCM
using JLD2
import DeltaRCM: save_checkpoint, load_checkpoint!

# Grids needed to continue a run exactly (padded arrays, halo included). Random numbers need no
# state: each parcel's stream derives from `(seed, timestep, parcel)`.
const RESTART_FIELDS = (:eta, :eta0, :eta_init, :depth, :stage, :uw, :ux, :uy, :qw, :qx, :qy,
                        :sand_frac, :active_layer)

function DeltaRCM.save_checkpoint(m, path::AbstractString)
    tmp = path * ".tmp"
    f = getfield(m, :f)
    jldopen(tmp, "w") do file
        file["time"] = m.time
        file["iter"] = m.iter
        file["H_SL"] = Float64(m.H_SL)
        file["seed"] = m.seed
        file["save_iter"] = m.save_iter
        file["since_data"] = m.since_data
        file["size"] = size(m)
        for name in RESTART_FIELDS
            file[String(name)] = getfield(f, name)
        end
    end
    mv(tmp, path; force = true)         # never leave a half-written checkpoint
    return nothing
end

function DeltaRCM.load_checkpoint!(m, path::AbstractString)
    f = getfield(m, :f)
    jldopen(path, "r") do file
        file["size"] == size(m) || throw(ArgumentError("Checkpoint grid $(file["size"]) does not match model $(size(m))"))
        m.time = file["time"]
        m.iter = file["iter"]
        m.H_SL = file["H_SL"]
        m.seed = file["seed"]
        m.save_iter = file["save_iter"]
        m.since_data = file["since_data"]
        m.since_checkpoint = 0.0
        for name in RESTART_FIELDS
            copyto!(getfield(f, name), file[String(name)])
        end
    end
    return nothing
end

end
