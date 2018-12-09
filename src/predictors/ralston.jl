export Ralston


"""
    Ralston()

The [Ralston](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Ralston's_method)
predictor of order 2.
"""
struct Ralston <: AbstractStatelessPredictor end
struct RalstonCache{T} <: AbstractStatelessPredictorCache
    A::Matrix{T}
    mk₁::Vector{T}
    mk₂::Vector{T}
end

function cache(::Ralston, H, x, t)
    mk₁ = dt(H, x, t)
    RalstonCache(jacobian(H, x, t), mk₁, copy(mk₁))
end
#
function predict!(xnext, ::Ralston, cache::RalstonCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    mk₁, mk₂ = cache.mk₁, cache.mk₂
    n = length(xnext)
    @inbounds for i=1:n
        mk₁[i] = -ẋ[i]
        xnext[i] = muladd(0.75 * Δt, ẋ[i], x[i])
    end

    minus_ẋ!(mk₂, H, xnext, t + 0.75 * Δt, cache.A)

    @inbounds for i=1:n
        xnext[i] = x[i] - Δt * (0.3333333333333333 * mk₁[i] + 0.6666666666666666 * mk₂[i])
    end

    nothing
end

order(::Ralston) = 3
