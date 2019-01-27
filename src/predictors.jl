export AbstractPredictor,
    AbstractPredictorCache,
    cache,
    setup!,
    update!,
    predict!

"""
    AbstractPredictor

The path tracking problem can be considered a initial value problem.
A predictor make use of this fact and predicts for a given pair ``(x₀,t)``
with ``H(x₀,t)=0`` a new ``x`` such that ``H(x, t + Δt)=0``.

The differential equation is
```math
Hₓ(x(t), t)ẋ(t)+H_t(x(t), t) = 0
```
which follows from the fact ``d/dt H(x(t),t) ≡ 0 ∀ t∈[0,1]`` and the total derivative
of ``H`` w.r.t. ``t``.
"""
abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache, x, ẋ, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
cache(p::AbstractPredictor, args...) = throw(MethodError(cache, tuple(p, args...)))


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t, Δt, ẋ)

Perform a prediction step for the value of `x` with step size `Δt`. `ẋ` is the derivate of `x(t)` at `t`.
"""
function predict! end

"""
    update!(cache::AbstractStatefulPredictorCache, H, x, ẋ, t, fac)

Update the cache. `x` is the new path value at `t` and `ẋ` is the derivative at `t`.
`fac` is a factorization of the Jacobian at `(x,t)`.
"""
update!(::AbstractPredictorCache, H, x, ẋ, t, fac) = nothing

"""
    setup!(cache::AbstractStatefulPredictorCache, H, x, ẋ, t, fac)

Setup the cache. `x` is the new path value at `t` and `ẋ` is the derivative at `t`.
`fac` is a factorization of the Jacobian at `(x,t)`. This falls back to calling `update`.
"""
setup!(C::AbstractPredictorCache, H, x, ẋ, t, fac) = update!(C, H, x, ẋ, t, fac)

"""
    reset!(cache, x, t)

Reset the cache. `x` is the path value at `t`.
"""
function reset!(::AbstractPredictorCache, x, t)
    nothing
end

"""
    order(::AbstractPredictor)

The order of a predictor is defined as the order of asymptotic error.
"""
function order end

# HELPERS

"""
    minus_ẋ!(out, H, x, t, J, dt)

Evaluate `Hₓ(x(t), t)⁻¹∂H∂t(x(t), t)` and store the result in `out`. `A` needs
to be able to store the Jacobian of `H`.
"""
function minus_ẋ!(ẋ, H, x, t, J, dt)
    jacobian_and_dt!(J, dt, H, x, t)
    solve!(ẋ, J, dt)
end

include("predictors/null_predictor.jl")
include("predictors/euler.jl")
include("predictors/rk3.jl")
include("predictors/rk4.jl")
include("predictors/heun.jl")
include("predictors/midpoint.jl")
include("predictors/ralston.jl")
include("predictors/pade21.jl")
