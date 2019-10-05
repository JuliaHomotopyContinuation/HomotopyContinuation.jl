# Core Tracker

A [`CoreTracker`](@ref) is the low-level path tracker. Its job is to track an initial solution
`x₁` from `t₁` to `t₀` by following a solution path ``x(t)``. For this a *predictor-corrector* scheme is used.
See our [introduction guide](https://www.juliahomotopycontinuation.org/guides/introduction/#tracking-solution-paths)
for a high-level explanation of the predictor corrector scheme


```@docs
CoreTracker
```

The easiest way to construct a `CoreTracker` is to use `coretracker` and `coretracker_startsolutions`.

```@docs
coretracker_startsolutions
coretracker
```

## Result and Status
```@docs
CoreTrackerResult
is_success(result::CoreTrackerResult)
solution(result::CoreTrackerResult)
CoreTrackerStatus.states
is_success(::CoreTrackerStatus.states)
is_terminated(::CoreTrackerStatus.states)
is_tracking(::CoreTrackerStatus.states)
is_invalid_startvalue(::CoreTrackerStatus.states)
```

## Methods
To track from a start to an endpoint with the `CoreTracker` we provide the following
routines.
```@docs
track(tracker::CoreTracker, x₁::AbstractVector, t₁::Number = 1.0, t₀::Number = 0.0)
track!(tracker::CoreTracker, x₁::AbstractVector, t₁::Number, t₀::Number)
init!(tracker::CoreTracker, x₁::AbstractVector, t₁::Number, t₀::Number)
```

It is also possible to use a `CoreTracker` as an iterator. This can either
be done by the high level `iterator` method or by directly using a `CoreTracker`
as an iterator. The recommend approach is to simply use `iterator`.
```@docs
iterator
```

## Introspecting the current state
To introspect the current state we provide the following routines.
```@docs
current_x
current_t
current_Δt
steps
status
LinearAlgebra.cond(::CoreTracker)
options
```

## Changing options
To change settings
```@docs
accuracy(::CoreTracker)
set_accuracy!
max_corrector_iters
set_max_corrector_iters!
max_step_size
set_max_step_size!
```

## Deprecated

The following functions are deprecated.
```@docs
max_refinement_iters
set_max_refinement_iters!
refinement_accuracy
set_refinement_accuracy!
```
