module Predictors


using ..Homotopies
using ..Utilities

export AbstractPredictor,
    AbstractPredictorCache,
    cache,
    predict!

"""
    AbstractPredictor

The path tracking problem can be considered a initial value problem.
A predictor make use of this fact and predicts for a given pair ``(x₀,t)``
with ``H(x₀,t)=0`` a new ``x`` such that ``H(x, t + Δt)=0``.

The differential equation is
```math
x′(t) = -Hₓ(x(t), t)⁻¹∂H∂t(x(t), t)
```
which follows from the fact ``∂tH(x(t),t) ≡ 0 ∀ t∈[0,1]`` and the total derivative
of ``H`` w.r.t. ``t``.
"""
abstract type AbstractPredictor end
abstract type AbstractStatelessPredictor <: AbstractPredictor end
abstract type AbstractStatefulPredictor <: AbstractPredictor end
abstract type AbstractPredictorCache end
abstract type AbstractStatelessPredictorCache <: AbstractPredictorCache end
abstract type AbstractStatefulPredictorCache <: AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache, x, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t, Δt)

Perform a prediction step for the value of `x` with step size `Δt`.
"""
function predict! end


setup!(::AbstractPredictorCache, x, ẋ, t) = nothing

"""
    update!(cache::AbstractStatefulPredictorCache, x, ẋ, t, Δt)

Update the cache. `x` is the new path value at `t+Δt` and `ẋ` is the derivative at `t + Δt`.
"""
update!(::AbstractPredictorCache, x, ẋ, t) = nothing

"""
    update_stepsize!(cache::AbstractStatefulPredictorCache, Δt_ratio)

Update the cache after a change of the step size.
`Δt_ratio` is the ratio between the new and the old step size.
"""
function update_stepsize!(::AbstractPredictorCache, Δt_ratio)
    nothing
end

"""
    reset!(cache, x, t)

Reset the cache. `x` is the path value at `t`.
"""
function reset!(::AbstractPredictorCache, x, t)
    nothing
end

function order end

function asymptotic_correction(alg::AbstractPredictor, ::AbstractPredictorCache,
                               correction_factor, Δs::Real)
    p = order(alg) + 1
    if p == 2
        √(correction_factor) * Δs
    else
        correction_factor^(1 // p) * Δs
    end
end

# HELPERS

"""
    minus_ẋ!(out, H, x, t, A)

Evaluate `Hₓ(x(t), t)⁻¹∂H∂t(x(t), t)` and store the result in `out`. `A` needs
to be able to store the Jacobian of `H`.
"""
function minus_ẋ!(out, H, x, t, A)
    jacobian_and_dt!(A, out, H, x, t)
    solve!(A, out)
end


include("predictors/null_predictor.jl")
include("predictors/euler.jl")
include("predictors/rk4.jl")
include("predictors/heun.jl")
include("predictors/midpoint.jl")
include("predictors/ralston.jl")

end
