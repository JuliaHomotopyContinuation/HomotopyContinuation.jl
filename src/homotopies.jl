module Homotopies

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import Base: start
import..Utilities: randomish_gamma


export AbstractHomotopy,
    StartTargetHomotopy,
    StraightLineHomotopy,
    γ,
    ParameterHomotopy,
    AbstractHomotopyCache,
    StraightLineHomotopyCache,
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

# StraightLineHomotopy implementation
include("homotopies/straight_line.jl")
include("homotopies/fixed_point.jl")
end
