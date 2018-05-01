module Homotopies

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
    nvariables,
    #start,
    target,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt

export HomotopyWithCache


"""
    AbstractHomotopy

Representing a homotopy.
"""
abstract type AbstractHomotopy end

"""
    AbstractStartTargetHomotopy

An homotopy with an explicit start and an explicit target system.
The systems should be of type [`AbstractSystem`](@ref).
"""
abstract type AbstractStartTargetHomotopy <: AbstractHomotopy end

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


abstract type AbstractParameterHomotopy <: AbstractHomotopy end


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
    jacobian_and_dt!(U, u, H, x, t [, cache::AbstractHomotopyCache])

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



"""
    nvariables(H::AbstractHomotopy)

Returns the number of variables of the homotopy `H`.
"""
nvariables(H::AbstractHomotopy) = last(size(H))


struct HomotopyWithCache{H<:AbstractHomotopy, C<:AbstractHomotopyCache} <: AbstractHomotopy
    homotopy::H
    cache::C
end

HomotopyWithCache(H::AbstractHomotopy, x, t) = HomotopyWithCache(H, cache(H, x, t))

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

function cache(H::AbstractStartTargetHomotopy, x, t)
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
struct StraightLineHomotopy{S<:AbstractSystem,T<:AbstractSystem} <: AbstractStartTargetHomotopy
    start::S
    target::T
    gamma::Complex{Float64}

end
function StraightLineHomotopy(start::AbstractSystem, target::AbstractSystem; gamma=randomish_gamma())
    StraightLineHomotopy(start, target, gamma)
end

function randomish_gamma()
    # Usually values near 1, i, -i, -1 are not good randomization
    # Therefore we artificially constrain the choices
    theta = rand() * 0.30 + 0.075 + (rand(Bool) ? 0.0 : 0.5)
    cis(2π * theta)
end

start(H::StraightLineHomotopy) = H.start
target(H::StraightLineHomotopy) = H.target

Base.size(H::StraightLineHomotopy) = size(H.start)

"""
    gamma(H::StraightLineHomotopy)

Obtain the gamma used in the StraightLineHomotopy.
"""
Base.Math.gamma(H::StraightLineHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the StraightLineHomotopy.
"""
γ(H::StraightLineHomotopy) = gamma(H)

function evaluate!(u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u

    u
end
function evaluate(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    (γ(H) * t) * G + (1 - t) * F
end
(H::StraightLineHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function dt!(u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= γ(H) .* c.u .- u

    u
end
function dt(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    γ(H) .* G .- F
end

function jacobian!(U, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.jacobian!(c.U, start(H), x, c.start)
    Systems.jacobian!(U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    U
end
function jacobian(H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    G = Systems.jacobian(start(H), x, c.start)
    F = Systems.jacobian(target(H), x, c.target)
    (γ(H) * t) .* G .+ (1 - t) .* F
end

function evaluate_and_jacobian!(u, U, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u
    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    nothing
end

function jacobian_and_dt!(U, u, H::StraightLineHomotopy, x, t, c::StartTargetHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U
    u .= γ(H) .* c.u .- u

    nothing
end

end
