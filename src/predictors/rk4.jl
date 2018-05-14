export RK4


"""
    RK4()

The classical [Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods)
predictor of order 4.
"""
struct RK4 <: AbstractPredictor end
struct RK4Cache{T} <: AbstractPredictorCache
    A::Matrix{T}
    k1::Vector{T}
    k2::Vector{T}
    k3::Vector{T}
    k4::Vector{T}
end

function cache(::RK4, H, x, t)
    k1 = dt(H, x, t)
    RK4Cache(jacobian(H, x, t), k1, copy(k1), copy(k1), copy(k1))
end
#
function predict!(xnext, ::RK4, cache::RK4Cache, H::HomotopyWithCache, x, t, Δt)
    mk₁, mk₂, mk₃, mk₄ = cache.k1, cache.k2, cache.k3, cache.k4

    minus_x_prime!(mk₁, H, x, t, cache.A)

    @. xnext = x - 0.5Δt * mk₁
    minus_x_prime!(mk₂, H, xnext, t + 0.5Δt, cache.A)

    @. xnext = x - 0.5Δt * mk₂
    minus_x_prime!(mk₃, H, xnext, t + 0.5Δt, cache.A)

    @. xnext = x - Δt * mk₃
    minus_x_prime!(mk₄, H, xnext, t + Δt, cache.A)

    @. xnext = x - 0.16666666666666666Δt * (mk₁ + 2mk₂ + 2mk₃ + mk₄)
    nothing
end
