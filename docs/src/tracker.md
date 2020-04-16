# Tracker

[`Tracker`](@ref) is a data structure to track for a given [`AbstractHomotopy`](@ref) ``H(x,t)`` a solution
``x`` from ``t₁ ∈ ℂ`` to ``t₀ ∈ ℂ``, i.e.,  ``H(x,t₁) = 0`` and ``x'`` with ``H(x',t₀) = 0`` is
returned.
This is done by following an implicitly defined smooth path ``x(t)`` using a predictor-corrector
scheme. In particular, it is assumed that for all ``t`` on the line segment between
``t₁`` and ``t₀`` the Jacobian ``H_x(x(t),t)`` has full column-rank.
The algorithm uses as an predictor a Padé approximant of order (2,1) and as a corrector
Newton's method. The details of the algorithm are described in the article [[Tim20]](https://arxiv.org/abs/1902.02968).

## Constructor and Options
```@docs
Tracker
TrackerOptions
TrackerParameters
DEFAULT_TRACKER_PARAMETERS
CONSERVATIVE_TRACKER_PARAMETERS
FAST_TRACKER_PARAMETERS
```

## Tracking

```@docs
track(tracker::Tracker, x₁::AbstractVector, t₁::Number, t₀::Number)
```

## Result
```@docs
TrackerResult
solution(::TrackerResult)
is_success(::TrackerResult)
is_invalid_startvalue(::TrackerResult)
steps(::TrackerResult)
accepted_steps(::TrackerResult)
rejected_steps(::TrackerResult)
```

## Low-level API


```@docs
track!(tracker::Tracker, x₁::AbstractVector, t₁::Number, t₀::Number)
init!(tracker::Tracker, ::TrackerResult, ::Number, ::Number)
TrackerCode
is_success(::TrackerCode.codes)
is_tracking(::TrackerCode.codes)
is_invalid_startvalue(::TrackerCode.codes)
is_terminated(::TrackerCode.codes)

```
