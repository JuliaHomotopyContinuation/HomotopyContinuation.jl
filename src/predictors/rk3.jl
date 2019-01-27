export RK3


"""
    RK3()

The classical [Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods)
predictor of order 3.
"""
struct RK3 <: AbstractPredictor end
struct RK3Cache{T} <: AbstractPredictorCache
    J::Matrix{T}
    dt::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    k3::Vector{T}
    k4::Vector{T}
end

function cache(::RK3, H, x, ẋ, t)
    RK3Cache(jacobian(H, x, t), dt(H, x, t), copy(ẋ), copy(ẋ), copy(ẋ), copy(ẋ))
end
#
function predict!(xnext, ::RK3, cache::RK3Cache, H::HomotopyWithCache, x, t, Δt, ẋ)
    J, dt, mk₁, mk₂, mk₃ = cache.J, cache.dt, cache.k1, cache.k2, cache.k3
    n = length(xnext)
    @inbounds for i=1:n
        mk₁[i] = -ẋ[i]
        xnext[i] = x[i] - 0.5Δt * mk₁[i]
    end
    minus_ẋ!(mk₂, H, xnext, t + 0.5Δt, J, dt)

    @inbounds for i=1:n
        xnext[i] = x[i] + Δt * mk₁[i] - 2Δt * mk₂[i]
    end
    minus_ẋ!(mk₃, H, xnext, t + Δt, J, dt)
    @inbounds for i=1:n
        xnext[i] = x[i] - Δt * mk₃[i]
    end

    @inbounds for i=1:n
        xnext[i] = x[i] - Δt/6 * (mk₁[i] + 4mk₂[i] + mk₃[i])
    end

    nothing
end

order(::RK3) = 4
