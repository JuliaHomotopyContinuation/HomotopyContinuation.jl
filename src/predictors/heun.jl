export Heun


"""
    Heun()

The [Heun](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Heun's_method)
predictor of order 2.
"""
struct Heun <: AbstractPredictor end
struct HeunCache{T} <: AbstractPredictorCache
    J::Matrix{T}
    dt::Vector{T}
    mk₁::Vector{T}
    mk₂::Vector{T}
end

function cache(::Heun, H, x, ẋ, t)
    HeunCache(jacobian(H, x, t), dt(H, x, t), copy(ẋ), copy(ẋ))
end
#
function predict!(xnext, ::Heun, cache::HeunCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    J, dt, mk₁, mk₂ = cache.J, cache.dt, cache.mk₁, cache.mk₂
    n = length(xnext)
    @inbounds for i=1:n
        mk₁[i] = -ẋ[i]
        xnext[i] = x[i] + Δt * ẋ[i]
    end

    minus_ẋ!(mk₂, H, xnext, t + Δt, J, dt)

    @inbounds for i=1:n
        xnext[i] = x[i] - 0.5 * Δt * (mk₁[i] + mk₂[i])
    end

    nothing
end

order(::Heun) = 3
