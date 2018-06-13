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
    n = length(xnext)
    minus_x_prime!(mk₁, H, x, t, cache.A)
    for i=1:n
        xnext[i] = x[i] - 0.5Δt * mk₁[i]
    end

    minus_x_prime!(mk₂, H, xnext, t + 0.5Δt, cache.A)
    for i=1:n
        xnext[i] = x[i] - 0.5Δt * mk₂[i]
    end

    minus_x_prime!(mk₃, H, xnext, t + 0.5Δt, cache.A)
    for i=1:n
        xnext[i] = x[i] - Δt * mk₃[i]
    end

    minus_x_prime!(mk₄, H, xnext, t + Δt, cache.A)
    for i=1:n
        xnext[i] = x[i] - 0.16666666666666666Δt * (mk₁[i] + 2mk₂[i] + 2mk₃[i] + mk₄[i])
    end

    nothing
end
