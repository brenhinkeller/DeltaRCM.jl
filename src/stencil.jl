# Cell types
const WALL = Int8(-2)    # land / impassable (inlet strip, plus the padding halo)
const EDGE = Int8(-1)    # open boundary: parcels leaving through it are exported
const OCEAN = Int8(0)    # ordinary interior cell (wet or dry)
const CHANNEL = Int8(1)  # inlet channel cells

# 3x3 neighbourhood in column-major order, so that direction `k` is element `k` of
# `A[i-1:i+1, j-1:j+1]`; `k == 5` is the cell itself. `i` is the downstream (first)
# index and `j` the cross-stream (second) index.
const DI = (-1, 0, 1, -1, 0, 1, -1, 0, 1)
const DJ = (-1, -1, -1, 0, 0, 0, 1, 1, 1)
const DIST = (√2, 1.0, √2, 1.0, 1.0, 1.0, √2, 1.0, √2)

# ---- random numbers -----------------------------------------------------------------------
# Every stochastic stream (one per parcel in water routing, one per parcel batch in
# sediment routing) is derived from `(seed, kind, timestep, ...)`, so a run is exactly
# reproducible, independent of how work is divided over threads, and a checkpoint needs
# only to store the seed and step counter.

@inline function splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"""Re-seed `rng` in place for the stream `(seed, a, b, c)` (allocation-free)."""
@inline function reseed!(rng::Random.Xoshiro, seed::UInt64, a::Integer, b::Integer, c::Integer)
    k = splitmix64(splitmix64(splitmix64(seed + UInt64(a)) + UInt64(b)) + UInt64(c))
    k2 = splitmix64(k)
    k3 = splitmix64(k2)
    k4 = splitmix64(k3)
    # Direct state assignment (Xoshiro256++ state words plus its fixed `s4` counter): a full
    # re-seed without allocating a new generator.
    rng.s0 = k; rng.s1 = k2; rng.s2 = k3; rng.s3 = k4; rng.s4 = 0x0000000000000000
    return rng
end
stream(seed::UInt64, a::Integer, b::Integer, c::Integer) = reseed!(Random.Xoshiro(), seed, a, b, c)

# stream "kinds"
const STREAM_WATER = 1
const STREAM_SAND = 2
const STREAM_MUD = 3
const STREAM_START = 4

"""
Pick index `k` with probability `w[k] / total`, where `total == sum(w)` is supplied.
Only strictly positive weights can be chosen, even under floating-point round-off.

Branch-free: the index is the number of cumulative sums that do not exceed the uniform draw
(plus one), which avoids the mispredicted early exit of a linear scan. Equivalent to
scanning the cumulative sums left to right.
"""
@inline function weighted_pick(rng, w::NTuple{9,T}, total::T) where {T}
    u = rand(rng, T) * total
    @inbounds begin
        c1 = w[1]; c2 = c1 + w[2]; c3 = c2 + w[3]; c4 = c3 + w[4]; c5 = c4 + w[5]
        c6 = c5 + w[6]; c7 = c6 + w[7]; c8 = c7 + w[8]; c9 = c8 + w[9]
        k = 1 + (u >= c1) + (u >= c2) + (u >= c3) + (u >= c4) + (u >= c5) + (u >= c6) + (u >= c7) + (u >= c8)
        k <= 9 && !(k == 9 && u >= c9) && return k
        # u reached the total through round-off: take the last positive weight
        for j in 9:-1:1
            w[j] > 0 && return j
        end
    end
    return 0
end

@inline wpow(x::T, θ) where {T} = θ == 1 ? x : θ == 2 ? x * x : x^T(θ)
