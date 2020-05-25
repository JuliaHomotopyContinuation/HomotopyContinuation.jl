export RK4


"""
    RK4()

The classical [Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods)
predictor of order 4.
"""
struct RK4 <: AbstractPredictor end
struct RK4Cache{T} <: AbstractPredictorCache
    dt::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    k3::Vector{T}
    k4::Vector{T}
end

function cache(::RK4, H, x, ẋ, t)
    RK4Cache(promote(dt(H, x, t), copy(ẋ), copy(ẋ), copy(ẋ), copy(ẋ))...)
end
#
function predict!(
    xnext,
    cache::RK4Cache,
    H::HomotopyWithCache,
    x,
    t,
    Δt,
    ẋ,
    Jac::JacobianMonitor,
)
    dt, mk₁, mk₂, mk₃, mk₄ = cache.dt, cache.k1, cache.k2, cache.k3, cache.k4
    n = length(xnext)
    @inbounds for i = 1:n
        mk₁[i] = -ẋ[i]
        xnext[i] = x[i] - 0.5Δt * mk₁[i]
    end

    minus_ẋ!(mk₂, H, xnext, t + 0.5Δt, Jac, dt)
    @inbounds for i = 1:n
        xnext[i] = x[i] - 0.5Δt * mk₂[i]
    end

    minus_ẋ!(mk₃, H, xnext, t + 0.5Δt, Jac, dt)
    @inbounds for i = 1:n
        xnext[i] = x[i] - Δt * mk₃[i]
    end

    minus_ẋ!(mk₄, H, xnext, t + Δt, Jac, dt)
    @inbounds for i = 1:n
        xnext[i] =
            x[i] - 0.16666666666666666Δt * (mk₁[i] + 2 * mk₂[i] + 2 * mk₃[i] + mk₄[i])
    end

    nothing
end

order(::RK4Cache) = 5
