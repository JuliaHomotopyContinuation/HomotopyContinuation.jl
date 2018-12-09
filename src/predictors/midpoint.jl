export Midpoint


"""
    Midpoint()

The (explicit) [Midpoint](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Explicit_midpoint_method)
predictor of order 2.
"""
struct Midpoint <: AbstractStatelessPredictor end
struct MidpointCache{T} <: AbstractStatelessPredictorCache
    A::Matrix{T}
    mk₂::Vector{T}
end

function cache(::Midpoint, H, x, t)
    mk₂ = dt(H, x, t)
    MidpointCache(jacobian(H, x, t), mk₂)
end
#
function predict!(xnext, ::Midpoint, cache::MidpointCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    mk₂ = cache.mk₂
    n = length(xnext)
    @inbounds for i=1:n
        xnext[i] = x[i] + 0.5Δt * ẋ[i]
    end

    minus_ẋ!(mk₂, H, xnext, t + 0.5Δt, cache.A)

    @inbounds for i=1:n
        xnext[i] = x[i] - Δt * mk₂[i]
    end

    nothing
end

order(::Midpoint) = 3
