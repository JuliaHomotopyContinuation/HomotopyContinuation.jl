export Euler

"""
    Euler()

This uses the explicit Euler method for prediction, also known as the
tangent predictor.
"""
struct Euler <: AbstractPredictor end
struct EulerCache{T} <: AbstractPredictorCache
    A::Matrix{T}
    b::Vector{T}
end

cache(::Euler, H, x, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))


function predict!(xnext, ::Euler, cache::EulerCache, H::HomotopyWithCache, x, t, Δt)
    minus_x_prime!(cache.b, H, x, t, cache.A)
    @. xnext = x - Δt * cache.b
    nothing
end
