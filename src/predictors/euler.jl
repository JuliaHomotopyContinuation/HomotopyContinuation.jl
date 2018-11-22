export Euler

"""
    Euler()

This uses the explicit Euler method for prediction, also known as the
tangent predictor.
"""
struct Euler <: AbstractStatelessPredictor end
struct EulerCache{T} <: AbstractStatelessPredictorCache
    A::Matrix{T}
    b::Vector{T}
end

cache(::Euler, H, x, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))

function predict!(xnext, ::Euler, cache::EulerCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    @inbounds for i=1:length(x)
        xnext[i] = x[i] + Δt * ẋ[i]
    end
    nothing
end

order(::Euler) = 1
