abstract type ParcelKind end
"Bedload parcels: transport-capacity limited deposition, erosion above `U_ero_sand`."
struct Sand <: ParcelKind end
"Suspended-load parcels: lag-limited deposition below `U_dep_mud`, erosion above `U_ero_mud`."
struct Mud <: ParcelKind end

"""
Running mass balance of sediment parcels (m³), reset each time step. Every parcel ends up
deposited on the bed, exported through an open boundary, or lost at `stepmax`.
"""
mutable struct SedimentTally{T}
    exported_sand::T
    exported_mud::T
    lost_sand::T
    lost_mud::T
end
SedimentTally{T}() where {T} = SedimentTally{T}(zero(T), zero(T), zero(T), zero(T))
function reset!(t::SedimentTally{T}) where {T}
    t.exported_sand = t.exported_mud = t.lost_sand = t.lost_mud = zero(T)
    return t
end

_theta(::Sand, c) = c.theta_sand
_theta(::Mud, c) = c.theta_mud

"""
    sediment_weights(kind, f, c, idx) -> (w, total)

Step probabilities for a sediment parcel: flow-direction weight times
`(depth * mod_sed_weight)^theta`, over wet non-wall neighbours. When all vanish the
parcel random-walks into any non-wall neighbour weighted by inverse distance.
"""
@inline function sediment_weights(kind::ParcelKind, f::Fields{T}, c::Constants{T}, idx::Int) where {T}
    off = c.off
    depth, ct = f.depth, f.cell_type
    qx = @inbounds f.qx[idx]
    qy = @inbounds f.qy[idx]
    dry = c.dry_depth
    θ = _theta(kind, c)
    w = ntuple(Val(9)) do k
        @inbounds n = idx + off[k]
        @inbounds if (k != 5) & (depth[n] > dry) & (ct[n] != WALL)
            max(zero(T), (qx * DI[k] + qy * DJ[k]) * T(INVDIST2[k])) * wpow(depth[n] * f.mod_sed_weight[n], θ)
        else
            zero(T)
        end
    end
    total = sum(w)
    total > 0 && return w, total
    w = ntuple(Val(9)) do k
        @inbounds (k != 5) & (ct[idx+off[k]] != WALL) ? T(INVDIST[k]) : zero(T)
    end
    return w, sum(w)
end

# erosion volume and the "no more than a quarter of the cell's water volume" limiter
@inline function limit_dVp(Vp::T, stage::T, eta::T, c::Constants{T}, mod_stable::T, erosion::Bool) where {T}
    d = stage - eta
    if d < 0
        erosion || return zero(T)
        d = -d
    end
    return min(Vp, mod_stable * d / 4 * c.dx^2)
end

@inline function erosion_volume(c::Constants{T}, U::T, U_ero::T, mod_erosion::T) where {T}
    return mod_erosion * c.Vp_sed * (U^c.beta - U_ero^c.beta) / U_ero^c.beta
end

"""
Apply a bed change `dVp` (m³, positive = deposition) at cell `n`, and refresh the local
depth and velocity.
"""
@inline function update_cell!(f::Fields{T}, c::Constants{T}, n::Int, dVp::T) where {T}
    @inbounds begin
        eta = f.eta[n] + dVp / c.dx^2
        d = max(f.stage[n] - eta, zero(T))
        f.eta[n] = eta
        f.depth[n] = d
        qw = f.qw[n]
        if d > c.dry_depth && qw > 0
            u = min(c.u_max, qw / d)
            f.uw[n] = u
            f.ux[n] = u * f.qx[n] / qw
            f.uy[n] = u * f.qy[n] / qw
        else
            f.uw[n] = 0; f.ux[n] = 0; f.uy[n] = 0
        end
    end
    return nothing
end

"""
    deposit_or_erode!(kind, f, c, n, Vp_res) -> Vp_res

Decide and apply deposition or erosion of a parcel carrying `Vp_res` at cell `n`;
returns the parcel's updated volume.
"""
@inline function deposit_or_erode!(::Sand, f::Fields{T}, c::Constants{T}, n::Int, Vp_res::T) where {T}
    @inbounds begin
        U = f.uw[n]
        qs_cap = c.qs0 * c.f_bedload / c.u0^c.beta * U^c.beta
        qs_loc = f.qs[n]
        dVp = zero(T)
        if qs_loc > qs_cap
            dVp = limit_dVp(Vp_res, f.stage[n], f.eta[n], c, f.mod_stable_weight[n], false)
            f.Vp_dep_sand[n] += dVp
        elseif U > c.U_ero_sand && qs_loc < qs_cap
            dVp = -limit_dVp(erosion_volume(c, U, c.U_ero_sand, f.mod_erosion[n]), f.stage[n], f.eta[n],
                             c, f.mod_stable_weight[n], true)
        end
        update_cell!(f, c, n, dVp)
        return Vp_res - dVp
    end
end

@inline function deposit_or_erode!(::Mud, f::Fields{T}, c::Constants{T}, n::Int, Vp_res::T) where {T}
    @inbounds begin
        U = f.uw[n]
        dVp = zero(T)
        if U < c.U_dep_mud
            dVp = c.lambda * Vp_res * (c.U_dep_mud^c.beta - U^c.beta) / c.U_dep_mud^c.beta
            dVp = limit_dVp(dVp, f.stage[n], f.eta[n], c, f.mod_stable_weight[n], false)
            f.Vp_dep_mud[n] += dVp
        elseif U > c.U_ero_mud
            dVp = -limit_dVp(erosion_volume(c, U, c.U_ero_mud, f.mod_erosion[n]), f.stage[n], f.eta[n],
                             c, f.mod_stable_weight[n], true)
        end
        update_cell!(f, c, n, dVp)
        return Vp_res - dVp
    end
end

# bedload carries a sediment-flux field `qs` that sets its local transport capacity check
@inline partition_flux!(::Mud, f, c, cur, nxt, Vp_res) = nothing
@inline function partition_flux!(::Sand, f::Fields{T}, c::Constants{T}, cur, nxt, Vp_res::T) where {T}
    part = Vp_res / 2 / c.dt / c.dx
    @inbounds f.qs[cur] += part
    @inbounds f.qs[nxt] += part
    return nothing
end

@inline _vp_dep(::Sand, f) = f.Vp_dep_sand
@inline _vp_dep(::Mud, f) = f.Vp_dep_mud
@inline _exported!(t::SedimentTally, ::Sand, v) = (t.exported_sand += v; nothing)
@inline _exported!(t::SedimentTally, ::Mud, v) = (t.exported_mud += v; nothing)
@inline _lost!(t::SedimentTally, ::Sand, v) = (t.lost_sand += v; nothing)
@inline _lost!(t::SedimentTally, ::Mud, v) = (t.lost_mud += v; nothing)

"""
    route_sediment_parcel!(kind, f, c, inlet, rng, tally)

Carry one sediment parcel from the inlet until it leaves through an open boundary or
reaches `stepmax`, depositing/eroding along the way. Sequential by nature: each
deposit changes the bed seen by the next parcel.
"""
function route_sediment_parcel!(kind::ParcelKind, f::Fields{T}, c::Constants{T}, inlet::Vector{Int},
                                rng, tally::SedimentTally{T}) where {T}
    (; Lp, off, stepmax) = c
    ct = f.cell_type
    Vp_res = c.Vp_sed
    cur = 2 + (inlet[rand(rng, 1:length(inlet))] - 1) * Lp
    kind isa Sand && (@inbounds f.qs[cur] += Vp_res / 2 / c.dt / c.dx)
    for it in 1:stepmax
        w, total = sediment_weights(kind, f, c, cur)
        if total <= 0
            _lost!(tally, kind, Vp_res)
            return nothing
        end
        nxt = cur + @inbounds off[weighted_pick(rng, w, total)]
        partition_flux!(kind, f, c, cur, nxt, Vp_res)
        cur = nxt
        if @inbounds ct[cur] == EDGE
            _exported!(tally, kind, Vp_res)
            return nothing
        end
        Vp_res = deposit_or_erode!(kind, f, c, cur, Vp_res)
        if it == stepmax
            if c.force_deposit
                @inbounds _vp_dep(kind, f)[cur] += Vp_res
                update_cell!(f, c, cur, Vp_res)
                Vp_res = zero(T)
            end
            _lost!(tally, kind, Vp_res)
        end
    end
    return nothing
end

"""
    route_sediment_parcels!(kind, f, c, inlet, seed, timestep, tally)

Route all parcels of `kind` for one time step (`Np_sed * f_bedload` sand parcels, the
rest mud).
"""
function route_sediment_parcels!(kind::ParcelKind, f::Fields{T}, c::Constants{T}, inlet::Vector{Int},
                                 seed::UInt64, timestep::Int, tally::SedimentTally{T}) where {T}
    nsand = floor(Int, c.Np_sed * c.f_bedload)
    n = kind isa Sand ? nsand : floor(Int, c.Np_sed * (1 - c.f_bedload))
    rng = stream(seed, kind isa Sand ? STREAM_SAND : STREAM_MUD, timestep, 0)
    for _ in 1:n
        route_sediment_parcel!(kind, f, c, inlet, rng, tally)
    end
    return nothing
end

"""
    topo_diffusion!(f, c)

Diffuse the bed between neighbouring open cells (walls excluded, so bed mass is conserved)
with a flux proportional to the summed sand flux of the pair and their elevation difference.
`N_crossdiff` Jacobi passes; the inlet row is held fixed.
"""
function topo_diffusion!(f::Fields{T}, c::Constants{T}) where {T}
    (; Lp, Wp, off) = c
    mult = c.diffusion_multiplier
    eta, qs, cf, mask, ct = f.eta, f.qs, f.cf, f.open_mask, f.cell_type
    for _ in 1:c.N_crossdiff
        @inbounds for j in 2:Wp-1, i in 2:Lp-1
            idx = i + (j - 1) * Lp
            if ct[idx] == WALL || i == 2
                cf[idx] = 0
            else
                e0, q0 = eta[idx], qs[idx]
                s = zero(T)
                for k in (1, 2, 3, 4, 6, 7, 8, 9)
                    (i == 3 && DI[k] == -1) && continue      # the inlet row is held fixed, so does not exchange
                    n = idx + off[k]
                    s += mask[n] * (q0 + qs[n]) * (eta[n] - e0)
                end
                cf[idx] = mult * s
            end
        end
        @inbounds for j in 2:Wp-1, i in 2:Lp-1
            eta[i, j] += cf[i, j]
        end
    end
    return nothing
end
