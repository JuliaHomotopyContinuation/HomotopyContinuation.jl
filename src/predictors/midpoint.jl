export Midpoint


"""
    Midpoint()

The (explicit) [Midpoint](https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods#Explicit_midpoint_method)
predictor of order 2.
"""
struct Midpoint <: AbstractPredictor end
struct MidpointCache{T} <: AbstractPredictorCache
    J::Matrix{T}
    dt::Vector{T}
    mk₂::Vector{T}
end

function cache(::Midpoint, H, x, ẋ, t)
    MidpointCache(jacobian(H, x, t), dt(H, x, t), copy(ẋ))
end
#
function predict!(xnext, ::Midpoint, cache::MidpointCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    J, dt, mk₂ = cache.J, cache.dt, cache.mk₂
    n = length(xnext)
    @inbounds for i=1:n
        xnext[i] = x[i] + 0.5Δt * ẋ[i]
    end

    minus_ẋ!(mk₂, H, xnext, t + 0.5Δt, J, dt)

    @inbounds for i=1:n
        xnext[i] = x[i] - Δt * mk₂[i]
    end

    nothing
end

order(::Midpoint) = 3
