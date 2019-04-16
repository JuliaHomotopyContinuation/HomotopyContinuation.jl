# Path Tracker

The [`solve`](@ref) routine is only a very thin wrapper around [`PathTracker`](@ref).
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


## Methods
To track a single path you can use the [`track`](@ref) and [`track!`](@ref) methods.
```@docs
track(tracker::PathTracker, x₁, t₁::Float64=1.0; path_number::Int=1, details::Symbol=:default, kwargs...)
track!(tracker::PathTracker, x₁, t₁::Float64=1.0; kwargs...)
```

The return type of [`track!`](@ref) is
```@docs
PathTrackerStatus.states
```

In the case that you track paths of parameter homotopy you can also change
the parameters using
```@docs
set_parameters!(::PathTracker)
```

## PathResult
For each path we return a [`PathResult`](@ref) containing the detailed information about
the single path.
```@docs
PathResult
```

The following helper functions are provided
```@docs
solution
residual
start_solution
Base.isreal(::PathResult)
LinearAlgebra.issuccess(::PathResult)
isfailed
isaffine
isprojective
isatinfinity
issingular
isnonsingular
```
