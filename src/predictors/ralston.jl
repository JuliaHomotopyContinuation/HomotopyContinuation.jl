export Ralston


"""
    Ralston()

The [Ralston](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Ralston's_method)
predictor of order 2.
"""
struct Ralston <: AbstractPredictor end
struct RalstonCache{T} <: AbstractPredictorCache
    dt::Vector{T}
    mk₁::Vector{T}
    mk₂::Vector{T}
end

function cache(::Ralston, H, x, ẋ, t)
    RalstonCache(dt(H, x,t), copy(ẋ), copy(ẋ))
end
#
function predict!(xnext, cache::RalstonCache, H::HomotopyWithCache, x, t, Δt, ẋ, Jac::JacobianMonitor)
    dt, mk₁, mk₂ = cache.dt, cache.mk₁, cache.mk₂
    n = length(xnext)
    @inbounds for i=1:n
        mk₁[i] = -ẋ[i]
        xnext[i] = muladd(0.75 * Δt, ẋ[i], x[i])
    end

    minus_ẋ!(mk₂, H, xnext, t + 0.75 * Δt, Jac, dt)

    @inbounds for i=1:n
        xnext[i] = x[i] - Δt * (0.3333333333333333 * mk₁[i] + 0.6666666666666666 * mk₂[i])
    end

    nothing
end

order(::RalstonCache) = 3
