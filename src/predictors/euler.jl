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

cache(::Euler, H, x, ẋ, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))

function predict!(xnext, cache::EulerCache, H::HomotopyWithCache, x, t, Δt, ẋ, ::JacobianMonitor)
    @inbounds for i=1:length(x)
        xnext[i] = x[i] + Δt * ẋ[i]
    end
    nothing
end

order(::EulerCache) = 2
