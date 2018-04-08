module Predictors

using ..NewHomotopies
using ..Utilities

export AbstractPredictor,
    AbstractPredictorCache,
    cache,
    predict!,
    NullPredictor, NullPredictorCache,
    Euler, EulerCache

abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache{M, N}, x, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t, dt)

Perform a prediction step for the value of `x` such that ``H(x, t+dt) ≈ 0``.
"""
function predict! end


struct NullPredictor <: AbstractPredictor end
struct NullPredictorCache <: AbstractPredictorCache end

cache(::NullPredictor, H, x, t) = NullPredictorCache()

function predict!(xnext, ::NullPredictor, ::NullPredictorCache, H, x, t, dt)
    xnext .= x
    nothing
end


struct Euler <: AbstractPredictor end
struct EulerCache{T} <: AbstractPredictorCache
    A::Matrix{T}
    b::Vector{T}
end

cache(::Euler, H, x, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))

function predict!(xnext, ::Euler, cache::EulerCache, H::HomotopyWithCache{N, N}, x, t, dt) where N
    A, b = cache.A, cache.b

    jacobian_and_dt!(A, b, H, x, t)
    solve_with_lu_inplace!(A, b)

    @. xnext = x - dt * b

    nothing
end


end
