# Path tracking

We also export a path tracking primitive to make the core path tracking routine
available for other applications.
At the heart is a [`PathTracker`](@ref) object which holds
all the state. The easiest way to construct a `PathTracker` is to use the [`pathtracker_startsolutions`](@ref) routine.

```@docs
pathtracker_startsolutions
pathtracker
```

## Types
```@docs
PathTracker
PathTrackerResult
PathTrackerStatus.states
```

## Methods
To track from a start to an endpoint with the `PathTracker` we provide the following
routines.
```@docs
track
track!
setup!
```

It is also possible to use a `PathTracker` as an iterator. This can either
be done by the high level `iterator` method or by directly using a `PathTracker`
as an iterator. The recommend approach is simply using `iterator`.
```@docs
iterator
```

## Introspecting the current state
To introspect the current state we provide the following routines.
```@docs
currx
currt
currΔt
curriters
currstatus
```

## Changing options
To change settings
```@docs
accuracy
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
