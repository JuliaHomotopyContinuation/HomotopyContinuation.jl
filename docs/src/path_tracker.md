# PathTracker

[`PathTracker`](@ref) is a data structure to track for a given [`AbstractHomotopy`](@ref) ``H(x,t)`` a solution
``x`` from ``t₁ > 0`` to ``0``, i.e.,  ``H(x,t₁) = 0`` and ``x'`` with ``H(x',0) = 0`` is
returned.
This is done by following an implicitly defined path ``x(t)`` using [`Tracker`](@ref).
In contrast to [`Tracker`](@ref) this uses an *endgame* to handle diverging paths and singular solutions.

```@docs
AbstractPathTracker
```

## Constructor and Options
```@docs
PathTracker
PathTrackerOptions
```


## Tracking

```@docs
track(::PathTracker, ::AbstractVector, ::Real)
```
