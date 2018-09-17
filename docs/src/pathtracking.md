# Path tracking

We also export a path tracking primitive to make the core path tracking routine
available for other applications.
At the heart is a [`PathTracking.PathTracker`](@ref) object which holds
all the state. The easiest way to construct a `PathTracker` is to use the [`pathtracker_startsolutions`](@ref) routine.

```@docs
pathtracker_startsolutions
```

## Types
```@docs
PathTracking.PathTracker
PathTracking.PathTrackerResult
StepLength.HeuristicStepLength
```

## Methods
To track from a start to an endpoint with the `PathTracker` we provide the following
routines.
```@docs
PathTracking.track!
PathTracking.track
PathTracking.setup!
```

To introspect the current state and change settings we provide the following routines.
```@docs
PathTracking.currx
PathTracking.currt
PathTracking.currΔt
PathTracking.curriters
PathTracking.currstatus
PathTracking.tol
PathTracking.corrector_maxiters
PathTracking.refinement_tol
PathTracking.refinement_maxiters
PathTracking.set_tol!
PathTracking.set_corrector_maxiters!
PathTracking.set_refinement_tol!
PathTracking.set_refinement_maxiters!
```
