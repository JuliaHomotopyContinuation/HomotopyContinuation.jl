module NewHomotopies

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import Base: start

export AbstractHomotopy,
    StartTargetHomotopy,
    StraightLineHomotopy,
    γ,
    ParameterHomotopy,
    AbstractHomotopyCache,
    StartTargetHomotopyCache,
    #start,
    target,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt

export HomotopyWithCache


"""
    AbstractHomotopy{M, N}

Representing a homotopy between `M` functions in `N` variables.
"""
abstract type AbstractHomotopy{M, N} end

"""
    AbstractStartTargetHomotopy

An homotopy with an explicit start and an explicit target system.
The systems should be of type [`AbstractSystem`](@ref).
"""
abstract type AbstractStartTargetHomotopy{M, N} <: AbstractHomotopy{M, N} end

"""
    start(::AbstractStartTargetHomotopy{M,N})::AbstractSystem{M,N}

Returns the start system.
"""
function start end

"""
    target(::AbstractStartTargetHomotopy{M,N})::AbstractSystem{M,N}

Returns the target system.
"""
function target end


abstract type AbstractParameterHomotopy{M, N} <: AbstractHomotopy{M, N} end


# Cache
"""
    AbstractHomotopyCache

A cache to avoid allocations for the evaluation of an [`AbstractHomotopy`](@ref).
"""
abstract type AbstractHomotopyCache end

"""
    cache(H::AbstractHomotopy, x, t)

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`.
"""
function cache end


# Homotopy API

"""
    evaluate!(u, H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
evaluate!(u, H::AbstractHomotopy, x, t, c=cache(H, x, t)) = evaluate!(u, H, x, t, c)

"""
    evaluate(H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` at `(x, t)`.
"""
evaluate(H::AbstractHomotopy, x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

"""
    dt!(u, H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
dt!(u, H::AbstractHomotopy, x, t, c=cache(H, x, t)) = dt!(u, H, x, t, c)

"""
    dt(H::AbstractHomotopy, x::AbstractVector [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` at `(x, t)`.
"""
dt(H::AbstractHomotopy, x, t, c=cache(H, x, t)) = dt(H, x, t, c)

"""
    jacobian!(u, H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the Jacobian of the homotopy `H` at `(x, t)` and store the result in `u`.
"""
jacobian!(U, H::AbstractHomotopy, x, t, c=cache(H, x, t)) = jacobian!(U, H, x, t, c)

"""
    jacobian(H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the Jacobian of the homotopy `H` at `(x, t)`.
"""
jacobian(H::AbstractHomotopy, x, t, c=cache(H, x, t)) = jacobian(H, x, t, c)


# Optional
"""
    evaluate_and_jacobian!(u, U, F, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` and its Jacobian at `(x, t)` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    evaluate!(u, H, x, t, c)
    jacobian!(U, H, x, t, c)
    nothing
end

"""
    evaluate_and_jacobian(H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` and its Jacobian at `(x, t)`.
"""
function evaluate_and_jacobian(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    u = evaluate(H, x, t, c)
    U = jacobian(H, x, t, c)
    u, U
end

"""
    evaluate_and_dt!(u, v, F, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)` and store the results in `u` (evalution)
and `v` (∂t).
"""
function jacobian_and_dt!(U, u, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    jacobian!(U, H, x, t, c)
    dt!(u, H, x, t, c)
    nothing
end

"""
    evaluate_and_dt(H::AbstractHomotopy, x, t [, cache::AbstractHomotopyCache])

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)`.
"""
function jacobian_and_dt(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    U = jacobian(H, x, t, c)
    u = dt(H, x, t, c)
    U, u
end


# derive
Base.size(H::AbstractHomotopy{M, N}) where {M, N} = (M, N)

struct HomotopyWithCache{M, N, H<:AbstractHomotopy{M, N}, C<:AbstractHomotopyCache} <: AbstractHomotopy{M, N}
    homotopy::H
    cache::C
end

# TODO: What to do with cache? Should we dispatch at the top on the NullCache similar to the Systems case?
function HomotopyWithCache(hom::H, cache::C) where {M, N, H<:AbstractHomotopy{M, N}, C<:AbstractHomotopyCache}
    HomotopyWithCache{M, N, H, C}(hom, cache)
end
HomotopyWithCache(hom::AbstractHomotopy, x, t) = HomotopyWithCache(hom, cache(H, x, t))

(H::HomotopyWithCache)(x, t) = evaluate(H, x, t)

evaluate(H::HomotopyWithCache, x, t) = evaluate(H.homotopy, x, t, H.cache)
evaluate!(u, H::HomotopyWithCache, x, t) = evaluate!(u, H.homotopy, x, t, H.cache)
jacobian(H::HomotopyWithCache, x, t) = jacobian(H.homotopy, x, t, H.cache)
jacobian!(U, H::HomotopyWithCache, x, t) = jacobian!(U, H.homotopy, x, t, H.cache)
dt(H::HomotopyWithCache, x, t) = dt(H.homotopy, x, t, H.cache)
dt!(u, H::HomotopyWithCache, x, t) = dt!(u, H.homotopy, x, t, H.cache)
jacobian_and_dt(H::HomotopyWithCache, x, t) = jacobian_and_dt(H.homotopy, x, t, H.cache)
jacobian_and_dt!(U, u, H::HomotopyWithCache, x, t) = jacobian_and_dt!(U, u, H.homotopy, x, t, H.cache)
evaluate_and_jacobian(H::HomotopyWithCache, x, t) = evaluate_and_jacobian(H.homotopy, x, t, H.cache)
evaluate_and_jacobian!(u, U, H::HomotopyWithCache, x, t) = evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.cache)


# Default caches
"""
    NullCache

An empty cache if no cache is necessary.
"""
struct NullCache <: AbstractHomotopyCache end

"""
    StartTargetHomotopyCache

An simple cache for `StartTargetHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct StartTargetHomotopyCache{SC, TC, T1, T2} <: AbstractHomotopyCache
    start::SC
    target::TC

    u::Vector{T1} # for evaluation
    U::Matrix{T2} # for Jacobian
end
const STHC{SC, TC, T1, T2} = StartTargetHomotopyCache{SC, TC, T1, T2}

function cache(H::AbstractStartTargetHomotopy{M, N}, x, t) where {M, N}
    start_cache = Systems.cache(H.start, x)
    target_cache = Systems.cache(H.target, x)

    u = Systems.evaluate(H.start, x)
    U = Systems.jacobian(H.start, x)

    StartTargetHomotopyCache(start_cache, target_cache, u, U)
end


# StraightLineHomotopy implementation
"""
    StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct StraightLineHomotopy{
    M, N,
    Start<:AbstractSystem{M, N},
    Target<:AbstractSystem{M, N}
    } <: AbstractStartTargetHomotopy{M, N}
    start::Start
    target::Target

    gamma::Complex{Float64}
end
function StraightLineHomotopy(start::S, target::T; gamma=cis(2π * rand())) where
    {M, N, S<:AbstractSystem{M, N}, T<:AbstractSystem{M, N}}

    StraightLineHomotopy{M, N, S, T}(start, target, gamma)
end

const SLH{M, N, Start, Target} = StraightLineHomotopy{M, N, Start, Target}

start(H::SLH) = H.start
target(H::SLH) = H.target

"""
    gamma(H::StraightLineHomotopy)

Obtain the gamma used in the StraightLineHomotopy.
"""
Base.Math.gamma(H::SLH) = H.gamma

"""
    γ(H)

Obtain the gamma used in the StraightLineHomotopy.
"""
γ(H::SLH) = gamma(H)

function evaluate!(u, H::SLH, x, t, c::STHC)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u

    u
end
function evaluate(H::SLH, x, t, c::STHC)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    (γ(H) * t) * G + (1 - t) * F
end

function dt!(u, H::SLH, x, t, c::STHC)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= γ(H) .* c.u .- u

    u
end
function dt(H::SLH, x, t, c::STHC)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    γ(H) .* G .- F
end

function jacobian!(U, H::SLH, x, t, c::STHC)
    Systems.jacobian!(c.U, start(H), x, c.start)
    Systems.jacobian!(U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    U
end
function jacobian(H::SLH, x, t, c::STHC)
    G = Systems.jacobian(start(H), x, c.start)
    F = Systems.jacobian(target(H), x, c.target)
    (γ(H) * t) .* G .+ (1 - t) .* F
end

function evaluate_and_jacobian!(u, U, H::SLH, x, t, c::STHC)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u
    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    nothing
end

function jacobian_and_dt!(U, u, H::SLH, x, t, c::STHC)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U
    u .= γ(H) .* c.u .- u

    nothing
end

end
