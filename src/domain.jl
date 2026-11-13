const GRAVITY = 9.81

"""
    Constants{T}

Immutable derived model constants (everything pyDeltaRCM computes in
`create_other_variables` / `create_boundary_conditions`), stored at precision `T`.
Grid indices `i` (downstream) and `j` (cross-stream) are 1-based; arrays carry a one-cell
halo, so the padded array size is `(L+2, W+2)` and cell `(i, j)` lives at `[i+1, j+1]`.
"""
struct Constants{T<:AbstractFloat}
    # grid
    L::Int
    W::Int
    Lp::Int                 # L + 2 (padded)
    Wp::Int                 # W + 2
    CTR::Int                # cross-stream centre column (1-based, unpadded)
    L0::Int                 # inlet length (cells)
    N0::Int                 # inlet width (cells)
    off::NTuple{9,Int}      # linear-index offsets of the 3x3 neighbourhood (padded array)
    dx::T
    # flow / sediment scales
    h0::T
    hb::T
    u0::T
    S0::T
    u_max::T
    C0::T
    dry_depth::T
    gamma::T
    V0::T
    Qw0::T
    qw0::T
    Qp_water::T
    qs0::T
    dVs::T
    Qs0::T
    Vp_sed::T
    dt::T
    # algorithm controls
    stepmax::Int
    itermax::Int
    Np_water::Int
    Np_sed::Int
    Nsmooth::Int
    Csmooth::T
    omega_sfc::T
    omega_flow::T
    omega_flow_iter::T
    theta_water::T
    theta_sand::T
    theta_mud::T
    beta::Int
    lambda::T
    f_bedload::T
    U_dep_mud::T
    U_ero_sand::T
    U_ero_mud::T
    N_crossdiff::Int
    diffusion_multiplier::T
    active_layer_thickness::T
    sand_frac_bc::T
    force_deposit::Bool
    SLR::T
    subsidence_rate::T
    start_subsidence::T
end

"""
    with_dt(c::Constants, dt) -> Constants

Constants for a time step of `dt` seconds: the sediment volume per step, the volume per
parcel and the diffusion multiplier all scale with `dt`.
"""
function with_dt(c::Constants{T}, dt) where {T}
    dt = T(dt)
    dVs = c.Qs0 * dt
    vals = ntuple(fieldcount(Constants)) do i
        n = fieldname(Constants, i)
        n === :dt ? dt :
        n === :dVs ? dVs :
        n === :Vp_sed ? dVs / c.Np_sed :
        n === :diffusion_multiplier ? c.diffusion_multiplier * dt / c.dt :
        getfield(c, i)
    end
    return Constants{T}(vals...)
end

"""
    resolve_parameters(p::DeltaParameters) -> DeltaParameters

Round `Length`, `Width`, `L0_meters` and `N0_meters` to whole numbers of cells (warning
if they change), exactly as pyDeltaRCM does.
"""
function resolve_parameters(p::DeltaParameters)
    dx = p.dx
    _round_to_cells(x) = isapprox(x / dx, round(x / dx); atol = 1e-9) ? x : round(x / dx) * dx
    Length = _round_to_cells(p.Length)
    Width = _round_to_cells(p.Width)
    Length == p.Length || @warn "Parameter `Length` was changed from $(p.Length) to $Length."
    Width == p.Width || @warn "Parameter `Width` was changed from $(p.Width) to $Width."
    L = round(Int, Length / dx)
    W = round(Int, Width / dx)
    L0 = max(1, min(round(Int, p.L0_meters / dx), L ÷ 4))
    N0 = max(3, min(round(Int, p.N0_meters / dx), W ÷ 4))
    L0_meters = L0 * dx
    N0_meters = N0 * dx
    L0_meters ≈ p.L0_meters || @warn "Parameter `L0_meters` was changed from $(p.L0_meters) to $L0_meters."
    N0_meters ≈ p.N0_meters || @warn "Parameter `N0_meters` was changed from $(p.N0_meters) to $N0_meters."
    return with(p; Length, Width, L0_meters, N0_meters)
end

function Constants{T}(p::DeltaParameters) where {T<:AbstractFloat}
    dx = p.dx
    L = round(Int, p.Length / dx)
    W = round(Int, p.Width / dx)
    Lp, Wp = L + 2, W + 2
    ctr0 = floor(Int, W / 2) - 1                      # pyDeltaRCM's 0-based centre
    ctr0 <= 1 && (ctr0 = floor(Int, W / 2))
    CTR = ctr0 + 1
    L0 = round(Int, p.L0_meters / dx)
    N0 = round(Int, p.N0_meters / dx)

    h0 = p.h0
    hb = something(p.hb, h0)
    u0 = p.u0
    C0 = p.C0_percent / 100
    dry_depth = min(0.1, 0.1 * h0)
    gamma = GRAVITY * p.S0 * dx / u0^2
    gamma > 0.1 && @warn "gamma = $(round(gamma; digits = 3)) exceeds the recommended maximum of 0.1; consider lowering S0 or dx, or raising u0."
    V0 = h0 * dx^2
    Qw0 = u0 * h0 * N0 * dx
    qw0 = u0 * h0
    Qp_water = Qw0 / p.Np_water
    qs0 = qw0 * C0
    dVs = 0.1 * N0^2 * V0
    Qs0 = Qw0 * C0
    Vp_sed = dVs / p.Np_sed
    stepmax = something(p.stepmax, 2 * (L + W))
    dt = dVs / Qs0
    N_crossdiff = max(1, round(Int, dVs / V0))
    diffusion_multiplier = dt / N_crossdiff * p.alpha * 0.5 / dx^2
    a = something(p.active_layer_thickness, h0 / 2)

    return Constants{T}(
        L, W, Lp, Wp, CTR, L0, N0, ntuple(k -> DI[k] + DJ[k] * Lp, Val(9)), dx,
        h0, hb, u0, p.S0, 2u0, C0, dry_depth, gamma, V0, Qw0, qw0, Qp_water, qs0, dVs, Qs0, Vp_sed, dt,
        stepmax, p.itermax, p.Np_water, p.Np_sed, p.Nsmooth, p.Csmooth, p.omega_sfc, p.omega_flow,
        2 / p.itermax, p.theta_water, p.coeff_theta_sand * p.theta_water, p.coeff_theta_mud * p.theta_water,
        p.beta, p.sed_lag, p.f_bedload, p.coeff_U_dep_mud * u0, p.coeff_U_ero_sand * u0,
        p.coeff_U_ero_mud * u0, N_crossdiff, diffusion_multiplier, a, p.sand_frac_bc,
        p.force_deposit, p.SLR, p.subsidence_rate, p.start_subsidence,
    )
end

"""
    Fields{T}

All gridded model state, as halo-padded `(L+2) × (W+2)` matrices. Arrays are only ever
mutated in place. `qx`/`ux` are the downstream (first-index) components and `qy`/`uy` the
cross-stream components.
"""
struct Fields{T<:AbstractFloat}
    cell_type::Matrix{Int8}
    open_mask::Matrix{T}        # 1 where a cell takes part in smoothing/diffusion (not WALL), 0 on walls and halo
    eta::Matrix{T}              # bed elevation
    eta0::Matrix{T}             # bed elevation at start of the current step
    eta_init::Matrix{T}
    stage::Matrix{T}
    depth::Matrix{T}
    qx::Matrix{T}
    qy::Matrix{T}
    qw::Matrix{T}
    ux::Matrix{T}
    uy::Matrix{T}
    uw::Matrix{T}
    qxn::Matrix{T}
    qyn::Matrix{T}
    qwn::Matrix{T}
    qs::Matrix{T}
    Vp_dep_sand::Matrix{T}
    Vp_dep_mud::Matrix{T}
    sand_frac::Matrix{T}
    active_layer::Matrix{T}
    sfc_visit::Matrix{T}
    sfc_sum::Matrix{T}
    Hnew::Matrix{T}
    Hsmth::Matrix{T}
    cf::Matrix{T}
    sigma::Matrix{T}            # subsidence per time step
    mod_water_weight::Matrix{T}
    mod_sed_weight::Matrix{T}
    mod_erosion::Matrix{T}
    mod_stable_weight::Matrix{T}
    sed_attract::Matrix{T}      # work array: per-cell (depth * mod_sed_weight)^theta for the sediment walk
end

const FIELD_NAMES = fieldnames(Fields)

const UNIT_FIELDS = (:mod_water_weight, :mod_sed_weight, :mod_erosion, :mod_stable_weight)

function Fields{T}(Lp::Int, Wp::Int) where {T}
    vals = map(FIELD_NAMES) do n
        n === :cell_type ? fill(WALL, Lp, Wp) :
        n in UNIT_FIELDS ? ones(T, Lp, Wp) : zeros(T, Lp, Wp)
    end
    return Fields{T}(vals...)
end

"""
Set the initial delta: inlet channel, land strip, ocean basin, free surface and flow
(`create_domain` in pyDeltaRCM). Returns the vector of padded inlet-cell column indices.
"""
function init_domain!(f::Fields{T}, c::Constants{T}) where {T}
    (; L, W, L0, N0, CTR, dx, h0, hb, S0, qw0) = c
    ctr0 = CTR - 1
    channel_start = ctr0 - round(Int, N0 / 2) + 1     # 0-based first channel column
    radius = min(L - 5, W / 2 - 5)

    for j in 1:W, i in 1:L
        land = i <= L0
        chan = land && channel_start < j <= channel_start + N0
        wet0 = !land                                   # ocean-type cell before edges are assigned
        isedge = i == L || j == 1 || j == W || hypot(i - 1 - 3, j - CTR) >= radius
        ct = land ? (chan ? CHANNEL : WALL) : (isedge ? EDGE : OCEAN)
        ip, jp = i + 1, j + 1
        f.cell_type[ip, jp] = ct
        f.open_mask[ip, jp] = ct == WALL ? zero(T) : one(T)
        f.stage[ip, jp] = land ? max(0, L0 - i) * dx * S0 : zero(T)
        f.depth[ip, jp] = chan ? h0 : wet0 ? hb : zero(T)
        f.qx[ip, jp] = chan ? qw0 : wet0 ? qw0 / 5 : zero(T)
    end
    @. f.qw = abs(f.qx)
    for idx in eachindex(f.depth)
        d = f.depth[idx]
        if d > 0
            f.ux[idx] = f.qx[idx] / d
            f.uy[idx] = f.qy[idx] / d
            f.uw[idx] = f.qw[idx] / d
        end
    end
    @. f.eta = f.stage - f.depth
    f.eta0 .= f.eta
    f.eta_init .= f.eta
    f.sand_frac .= c.sand_frac_bc
    f.active_layer .= c.sand_frac_bc
    return [j + 1 for j in 1:W if f.cell_type[2, j+1] == CHANNEL]
end

"""Scratch space for [`flooding_correction!`](@ref) (reused to avoid allocation)."""
struct FloodBuffer{T}
    idx::Vector{Int}
    val::Vector{T}
end
FloodBuffer{T}() where {T} = FloodBuffer{T}(Int[], T[])

