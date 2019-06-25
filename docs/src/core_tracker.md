# Core Tracker

We also export the path tracking primitive to make the core path tracking routine
available for other applications.
At the heart is a [`CoreTracker`](@ref) object which holds
all the state.

```@docs
CoreTracker
```

The easiest way to construct a `CoreTracker`:
```@docs
coretracker_startsolutions
coretracker
```

## Other types
```@docs
CoreTrackerResult
CoreTrackerStatus.states
CoreTrackerOptions
```

## Methods
To track from a start to an endpoint with the `CoreTracker` we provide the following
routines.
```@docs
track(tracker::CoreTracker, x₁::AbstractVector, t₁=1.0, t₀=0.0; kwargs...)
track!(x₀, tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0; setup_patch=tracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)
setup!
```

It is also possible to use a `CoreTracker` as an iterator. This can either
be done by the high level `iterator` method or by directly using a `CoreTracker`
as an iterator. The recommend approach is simply using `iterator`.
```@docs
iterator
```

## Introspecting the current state
To introspect the current state we provide the following routines.
```@docs
current_x
current_t
current_Δt
iters
status
LinearAlgebra.cond(::CoreTracker)
digits_lost
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
max_refinement_iters
set_max_refinement_iters!
refinement_accuracy
set_refinement_accuracy!
```
