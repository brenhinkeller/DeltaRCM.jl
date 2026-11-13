const INVDIST = map(inv, DIST)
const INVDIST2 = map(d -> inv(d^2), DIST)
const UI = map(/, DI, DIST)     # downstream component of the unit direction vector
const UJ = map(/, DJ, DIST)     # cross-stream component

# ---- water-routing work buffers -----------------------------------------------------------

"""
Work space of the water parcel walk: it aliases the model's own accumulation arrays
`qxn, qyn, qwn, sfc_sum, sfc_visit`, and adds the loop-detection stamps and path buffers.
(A threaded variant that kept one of these per chunk is described in `notes/`.)
"""
struct WaterScratch{T}
    qxn::Matrix{T}
    qyn::Matrix{T}
    qwn::Matrix{T}
    sfc_sum::Matrix{T}
    sfc_visit::Matrix{T}
    stamp::Matrix{Int32}        # loop detection: `stamp[cell] == id` iff this parcel has been there
    path_i::Vector{Int32}
    path_j::Vector{Int32}
    counter::Base.RefValue{Int32}
end

function WaterScratch(f::Fields{T}, c::Constants{T}) where {T}
    return WaterScratch{T}(
        f.qxn, f.qyn, f.qwn, f.sfc_sum, f.sfc_visit,
        zeros(Int32, c.Lp, c.Wp), zeros(Int32, c.stepmax + 2), zeros(Int32, c.stepmax + 2), Ref(Int32(0)),
    )
end

function reset!(ch::WaterScratch)
    fill!(ch.qxn, 0); fill!(ch.qyn, 0); fill!(ch.qwn, 0)
    fill!(ch.sfc_sum, 0); fill!(ch.sfc_visit, 0)
    return ch
end

# ---- neighbour weights --------------------------------------------------------------------

"""
    water_weights(f, c, idx) -> (w, total)

Unnormalised probabilities of a water parcel at (padded, linear) cell `idx` stepping to
each of its 9 neighbours (`w[5] == 0`): a `gamma`-blend of the surface-slope and
flow-direction weights (each normalised), scaled by `(depth * mod)^theta_water`. Dry cells
and walls are never entered. If every weight vanishes the parcel random-walks into any
wet neighbour, then into any non-wall neighbour; `total == 0` means it is boxed in.
"""
@inline function water_weights(f::Fields{T}, c::Constants{T}, idx::Int) where {T}
    off = c.off
    stage, depth, ct = f.stage, f.depth, f.cell_type
    s0 = @inbounds stage[idx]
    qx = @inbounds f.qx[idx]
    qy = @inbounds f.qy[idx]
    dry = c.dry_depth
    ok = ntuple(Val(9)) do k
        @inbounds n = idx + off[k]
        @inbounds (k != 5) & (depth[n] > dry) & (ct[n] != WALL)
    end
    ws = ntuple(Val(9)) do k
        @inbounds ok[k] ? max(zero(T), (s0 - stage[idx+off[k]]) * T(INVDIST[k])) : zero(T)
    end
    wi = ntuple(Val(9)) do k
        ok[k] ? max(zero(T), (qx * DI[k] + qy * DJ[k]) * T(INVDIST2[k])) : zero(T)
    end
    ss, si = sum(ws), sum(wi)
    ns = ss > 0 ? inv(ss) : zero(T)
    ni = si > 0 ? inv(si) : zero(T)
    γ = c.gamma
    θ = c.theta_water
    w = ntuple(Val(9)) do k
        if ok[k]
            @inbounds n = idx + off[k]
            @inbounds max(zero(T), (γ * ws[k] * ns + (1 - γ) * wi[k] * ni) * wpow(depth[n] * f.mod_water_weight[n], θ))
        else
            zero(T)
        end
    end
    total = sum(w)
    total > 0 && return w, total

    nok = count(ok)
    if nok > 0
        return map(o -> o ? one(T) : zero(T), ok), T(nok)
    end
    open = ntuple(k -> @inbounds((k != 5) & (ct[idx+off[k]] != WALL)), Val(9))
    nopen = count(open)
    return map(o -> o ? one(T) : zero(T), open), T(nopen)
end

# ---- one water parcel ---------------------------------------------------------------------

"""
Route one water parcel from the inlet to an open boundary (or `stepmax`), accumulating its
flux into `qxn, qyn, qwn`, and -- if it reached the boundary without being
disqualified by a loop -- its contribution to the free surface into `sfc_sum, sfc_visit`
(the parcel's path is walked backwards from the boundary, where stage equals `H_SL`).
"""
function route_water_parcel!(ch::WaterScratch{T}, f::Fields{T}, c::Constants{T}, inlet::Vector{Int},
                             rng, H_SL::T) where {T}
    (; Lp, L, W, L0, CTR, off, stepmax, dx) = c
    ct = f.cell_type
    half_q = c.Qp_water / dx / 2

    # id for loop detection (never 0, resets on overflow)
    id = ch.counter[] + Int32(1)
    if id == typemax(Int32)
        fill!(ch.stamp, 0)
        id = Int32(1)
    end
    ch.counter[] = id

    # ghost-node flux into the inlet cell
    j = inlet[rand(rng, 1:length(inlet))]
    ci = 2
    cj = j
    cur = ci + (cj - 1) * Lp
    @inbounds begin
        ch.qxn[cur] += 1
        ch.qwn[cur] += half_q
        ch.path_i[1] = ci
        ch.path_j[1] = cj
    end
    n = 1
    flagged = false          # looped away from sea level: excluded from the free surface
    reached = false

    @inbounds while n <= stepmax
        w, total = water_weights(f, c, cur)
        total > 0 || break
        k = weighted_pick(rng, w, total)
        ni = ci + DI[k]
        nj = cj + DJ[k]
        nxt = cur + off[k]

        ux, uy = T(UI[k]), T(UJ[k])
        ch.qxn[cur] += ux; ch.qxn[nxt] += ux
        ch.qyn[cur] += uy; ch.qyn[nxt] += uy
        ch.qwn[cur] += half_q; ch.qwn[nxt] += half_q

        # loops: revisiting a cell below the inlet strip relocates the parcel outward
        if ni > L0 + 1
            if ch.stamp[nxt] == id
                abs(f.stage[nxt] - H_SL) < T(0.1) || (flagged = true)
                px0, py0 = ni - 2, nj - 2                # 0-based indices
                Fx, Fy = px0 - 1, py0 - (CTR - 1)
                Fw = hypot(Fx, Fy)
                px, py = px0, py0
                if Fw != 0
                    px += round(Int, Fx / Fw * 5)
                    py += round(Int, Fy / Fw * 5)
                end
                px = min(L - 2, max(px, L0))
                py = min(W - 2, max(1, py))
                ni, nj = px + 2, py + 2
                nxt = ni + (nj - 1) * Lp
            end
            ch.stamp[nxt] = id
        end

        n += 1
        ch.path_i[n] = ni
        ch.path_j[n] = nj
        cur, ci, cj = nxt, ni, nj
        if ct[cur] == EDGE
            reached = true
            break
        end
    end

    (reached && !flagged) || return nothing

    # free surface: walk back from the boundary (stage = H_SL) up the parcel's path
    u0, h0, S0 = c.u0, c.h0, c.S0
    H = H_SL
    dH = zero(T)
    in_ocean = true
    @inbounds for t in (n-1):-1:1
        i, j = Int(ch.path_i[t]), Int(ch.path_j[t])
        ip, jp = Int(ch.path_i[t+1]), Int(ch.path_j[t+1])
        idx = i + (j - 1) * Lp
        if i != ip || j != jp
            if in_ocean
                if f.uw[idx] > u0 / 2 || f.depth[idx] < h0 / 10
                    in_ocean = false
                end
            else
                uw = f.uw[idx]
                if uw == 0 || abs(ip - i) > 1 || abs(jp - j) > 1
                    dH = zero(T)
                else
                    dH = S0 * (f.ux[idx] * (ip - i) * dx + f.uy[idx] * (jp - j) * dx) / uw
                end
            end
        end
        H += dH
        ch.sfc_visit[idx] += 1
        ch.sfc_sum[idx] += H
    end
    return nothing
end

# ---- one water iteration ------------------------------------------------------------------

"""
    run_water_parcels!(scratch, f, c, inlet, seed, timestep, iteration, H_SL)

Route all `Np_water` parcels of one water iteration. Parcels do not change the bed or flow
they read, only accumulate into the flux and free-surface sums. Each parcel draws from its own
random stream, derived from `(seed, timestep, iteration, parcel)`.
"""
function run_water_parcels!(scratch::WaterScratch{T}, f::Fields{T}, c::Constants{T}, inlet::Vector{Int},
                            seed::UInt64, timestep::Int, iteration::Int, H_SL::T) where {T}
    Np = c.Np_water
    reset!(scratch)
    rng = Random.Xoshiro()
    for p in 1:Np
        reseed!(rng, seed, STREAM_WATER, timestep, (iteration - 1) * Np + p)
        route_water_parcel!(scratch, f, c, inlet, rng, H_SL)
    end
    return nothing
end

"""
Combine the parcel-accumulated surface into a new, smoothed, under-relaxed `stage`.
"""
function compute_free_surface!(f::Fields{T}, c::Constants{T}, H_SL::T, fb::FloodBuffer{T}) where {T}
    (; Lp, Wp) = c
    Hnew = f.Hnew
    @inbounds for j in 2:Wp-1, i in 2:Lp-1
        v = f.sfc_visit[i, j]
        h = v > 0 ? f.sfc_sum[i, j] / v : f.eta[i, j] + f.depth[i, j]
        Hnew[i, j] = max(h, H_SL, f.eta[i, j])
    end
    Hs = smooth_free_surface!(f, c)
    ω = c.omega_sfc
    @inbounds for j in 2:Wp-1, i in 2:Lp-1
        f.stage[i, j] = (1 - ω) * f.stage[i, j] + ω * Hs[i, j]
    end
    flooding_correction!(f, c, fb)
    return nothing
end

"""
Smooth `f.Hnew` with `Nsmooth` Jacobi passes of the centre-weighted neighbour average
(walls excluded, boundary cells held fixed). Returns whichever of `Hnew`/`Hsmth`
holds the result.
"""
function smooth_free_surface!(f::Fields{T}, c::Constants{T}) where {T}
    (; Lp, Wp, off, Nsmooth) = c
    src, dst = f.Hnew, f.Hsmth
    copyto!(dst, src)
    C = c.Csmooth
    ct, mask = f.cell_type, f.open_mask
    for _ in 1:Nsmooth
        @inbounds for j in 2:Wp-1, i in 2:Lp-1
            idx = i + (j - 1) * Lp
            t = ct[idx]
            if t != WALL && t != EDGE
                s = zero(T); cnt = zero(T)
                for k in (1, 2, 3, 4, 6, 7, 8, 9)
                    n = idx + off[k]
                    m = mask[n]
                    s += m * src[n]
                    cnt += m
                end
                dst[idx] = cnt > 0 ? C * src[idx] + (1 - C) * s / cnt : src[idx]
            end
        end
        src, dst = dst, src
    end
    return src
end

"""
Flood dry cells whose wet neighbours stand above the local bed: the cell takes the highest
wet-neighbour stage. Wetness is judged from `stage - eta` and all updates use the
pre-pass state (so the result does not depend on loop order).
"""
function flooding_correction!(f::Fields{T}, c::Constants{T}, fb::FloodBuffer{T}) where {T}
    (; Lp, Wp, off, dry_depth) = c
    stage, eta, ct = f.stage, f.eta, f.cell_type
    empty!(fb.idx); empty!(fb.val)
    @inbounds for j in 2:Wp-1, i in 2:Lp-1
        idx = i + (j - 1) * Lp
        (ct[idx] == WALL || stage[idx] - eta[idx] > dry_depth) && continue
        top = typemin(T)
        for k in (1, 2, 3, 4, 6, 7, 8, 9)
            n = idx + off[k]
            if ct[n] != WALL && stage[n] - eta[n] > dry_depth
                top = max(top, stage[n])
            end
        end
        if top > eta[idx]
            push!(fb.idx, idx); push!(fb.val, top)
        end
    end
    @inbounds for (idx, v) in zip(fb.idx, fb.val)
        stage[idx] = v
    end
    return nothing
end

"""
Recompute depth from stage, then the flow field (parcel flux direction scaled to the
accumulated discharge, under-relaxed) and velocities.
"""
function finalize_water_iteration!(f::Fields{T}, c::Constants{T}, inlet::Vector{Int},
                                   iteration::Int, timestep::Int, H_SL::T) where {T}
    (; Lp, Wp, dry_depth, u_max, qw0) = c
    ω = iteration == 1 ? c.omega_flow : c.omega_flow_iter
    first_step = timestep == 0
    @inbounds for j in 2:Wp-1, i in 2:Lp-1
        s = max(f.stage[i, j], H_SL)
        f.stage[i, j] = s
        d = max(s - f.eta[i, j], zero(T))
        f.depth[i, j] = d

        qxn, qyn, qwn = f.qxn[i, j], f.qyn[i, j], f.qwn[i, j]
        dloc = hypot(qxn, qyn)
        if dloc > 0
            r = qwn / dloc
            qxn *= r; qyn *= r
        end
        if first_step
            qx, qy = qxn, qyn
        else
            qx = ω * qxn + (1 - ω) * f.qx[i, j]
            qy = ω * qyn + (1 - ω) * f.qy[i, j]
        end
        f.qx[i, j] = qx
        f.qy[i, j] = qy
        f.qw[i, j] = hypot(qx, qy)
    end
    @inbounds for jp in inlet
        f.qx[2, jp] = qw0
        f.qy[2, jp] = 0
        f.qw[2, jp] = qw0
    end
    @inbounds for j in 2:Wp-1, i in 2:Lp-1
        d, qw = f.depth[i, j], f.qw[i, j]
        if d > dry_depth && qw > 0
            u = min(u_max, qw / d)
            f.uw[i, j] = u
            f.ux[i, j] = u * f.qx[i, j] / qw
            f.uy[i, j] = u * f.qy[i, j] / qw
        else
            f.uw[i, j] = 0; f.ux[i, j] = 0; f.uy[i, j] = 0
        end
    end
    return nothing
end
