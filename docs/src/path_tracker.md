# Path Tracker

The [`solve`](@ref) routine is only a thin wrapper around [`PathTracker`](@ref).
Therefore you can also use [`PathTracker`](@ref) directly.
This is for example a good choice if you have to solve the same problem many times.

```@docs
PathTracker
```

The easiest way to construct a `PathTracker`:
```@docs
pathtracker_startsolutions
pathtracker
```
To track a single path you can use the [`track`](@ref) and [`track!`](@ref) methods.
```@docs
track(tracker::PathTracker, x₁)
```

## PathResult
For each path we return a [`PathResult`](@ref) containing the detailed information about
the single path.
```@docs
PathResult
```

The following helper functions are provided
```@docs
solution(::PathResult)
accuracy(::PathResult)
residual(::PathResult)
winding_number(tracker::PathTracker)
multiplicity(::PathResult)
condition_jacobian(::PathResult)
LinearAlgebra.cond(::PathResult)
start_solution(::PathResult)
is_real(::PathResult)
is_success(::PathResult)
is_failed(::PathResult)
is_affine(::PathResult)
is_projective(r::PathResult)
is_at_infinity(::PathResult)
is_singular(::PathResult)
is_nonsingular(::PathResult)
```

## Low-level API

```@docs
track!(tracker::PathTracker, x₁)
```

In the case that you track paths of parameter homotopy you can also change
the parameters using
```@docs
start_parameters!(::PathTracker, p)
target_parameters!(::PathTracker, p)
```

The return type of [`track!`](@ref) is a [`PathTrackerStatus.states`](@ref):
```@docs
PathTrackerStatus.states
is_success
is_at_infinity
is_invalid_startvalue
is_failed
is_terminated_callback
is_tracking
```
