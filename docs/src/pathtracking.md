# Path tracking

We also export a path tracking primitive to make the core path tracking routine
available for other applications.
At the heart is a [`PathTracking.PathTracker`](@ref) object which holds
all the state. The easiest way to construct a `PathTracker` is to use the [`pathtracker_startsolutions`](@ref) routine.

```@docs
pathtracker_startsolutions
pathtracker
```

## Types
```@docs
PathTracker
PathTrackerResult
```

## Methods
To track from a start to an endpoint with the `PathTracker` we provide the following
routines.
```@docs
PathTracking.track
PathTracking.track!
PathTracking.setup!
```

It is also possible to use a `PathTracker` as an iterator:
```
iterator!
```

To introspect the current state and change settings we provide the following routines.
```@docs
currx
currt
currΔt
curriters
currstatus
tol
corrector_maxiters
refinement_tol
refinement_maxiters
set_tol!
set_corrector_maxiters!
set_refinement_tol!
set_refinement_maxiters!
```
