module Predictors

using ..NewHomotopies
using ..Utilities

export AbstractPredictor,
    AbstractPredictorCache,
    cache,
    predict!,
    NullPredictor, NullPredictorCache,
    Euler, EulerCache

"""
    AbstractPredictor

The path tracking problem can be considered a initial value problem.
A predictor make use of this fact and predicts for a given pair ``(x₀,t)``
with ``H(x₀,t)=0`` a new ``x`` such that ``H(x, t + Δt)=0``.

The differential equation is
```math
x′(t) = -Hₓ(x(t), t)⁻¹Hₜ(x(t), t)
```
which follows from the fact ``∂tH(x(t),t) ≡ 0 ∀ t∈[0,1]`` and the total derivative
of ``H`` w.r.t. ``t``.
"""
abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache{M, N}, x, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t, Δt)

Perform a prediction step for the value of `x` with step size `Δt`.
"""
function predict! end


struct NullPredictor <: AbstractPredictor end
struct NullPredictorCache <: AbstractPredictorCache end

cache(::NullPredictor, H, x, t) = NullPredictorCache()

function predict!(xnext, ::NullPredictor, ::NullPredictorCache, H, x, t, dt)
    xnext .= x
    nothing
end


struct Euler <: AbstractPredictor end
struct EulerCache{T} <: AbstractPredictorCache
    A::Matrix{T}
    b::Vector{T}
end

cache(::Euler, H, x, t) = EulerCache(jacobian(H, x, t), dt(H, x, t))

"""
    minus_x_prime!(out, H, x, t, A)

Evaluate `Hₓ(x(t), t)⁻¹Hₜ(x(t), t)` and store the result in `out`. `A` needs
to be able to store the Jacobian of `H`.
"""
function minus_x_prime!(out, H, x, t, A)
    jacobian_and_dt!(A, out, H, x, t)
    solve_with_lu_inplace!(A, out)
    out
end


function predict!(xnext, ::Euler, cache::EulerCache, H::HomotopyWithCache{N, N}, x, t, Δt) where N
    minus_x_prime!(cache.b, H, x, t, cache.A)
    @. xnext = x - Δt * cache.b
    nothing
end

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
    RK4Cache(jacobian(H, x, t), k1, copy(k1), copy(k1), copy(k1), copy(k1))
end
#
function predict!(xnext, ::RK4, cache::RK4Cache, H::HomotopyWithCache{N, N}, x, t, Δt) where N
    mk₁, mk₂, mk₃, mk₄ = cache.k1, cache.k2, cache.k3, cache.k4

    minus_x_prime!(mk₁, H, x, t, cache.A)

    @. xnext = x - 0.5Δt * mk₁
    minus_x_prime!(mk₂, H, xnext, t + 0.5Δt, cache.A)

    @. xnext .= x - 0.5Δt * mk₂
    minus_x_prime!(mk₃, H, xnext, t + 0.5Δt, cache.A)

    @. xnext = x - Δt * mk₃
    minus_x_prime!(mk₄, H, xnext, t + Δt, cache.A)

    @. xnext = x - 0.16666666666666666Δt * (mk₁ + 2mk₂ + 2mk₃ + mk₄)
    nothing
end
end
