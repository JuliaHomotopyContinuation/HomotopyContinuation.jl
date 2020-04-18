# PathTracker

[`PathTracker`](@ref) is a data structure to track for a given [`AbstractHomotopy`](@ref) ``H(x,t)`` a solution
``x`` from ``t₁ > 0`` to ``0``, i.e.,  ``H(x,t₁) = 0`` and ``x'`` with ``H(x',0) = 0`` is
returned.
This is done by following an implicitly defined path ``x(t)`` using [`Tracker`](@ref).
In contrast to [`Tracker`](@ref) this uses an *endgame* to handle diverging paths and singular solutions.


## Constructor and Options
```@docs
PathTracker
PathTrackerOptions
```


## Tracking

```@docs
track(::PathTracker, ::AbstractVector, ::Real)
```

## PathResult

```@docs
PathResult
solution(::PathResult)
is_success(::PathResult)
is_at_infinity(::PathResult)
is_excess_solution(::PathResult)
is_failed(::PathResult)
is_finite(::PathResult)
is_singular(::PathResult)
is_nonsingular(::PathResult)
is_real(::PathResult)
accuracy(::PathResult)
residual(::PathResult)
steps(::PathResult)
accepted_steps(::PathResult)
rejected_steps(::PathResult)
winding_number(::PathResult)
path_number(::PathResult)
start_solution(::PathResult)
multiplicity(::PathResult)
last_path_point(::PathResult)
valuation(::PathResult)
```
