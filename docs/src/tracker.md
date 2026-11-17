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

## Extended precision

If the option `extended_precision` of [`TrackerOptions`](@ref) is `true` (the default), the tracker switches to double-double arithmetic (about 32 digits), if iterative refinement with Float64 precision fails. 
The path coordinates, `t` and the returned solutions always remain `ComplexF64`. If the option is `false`, only `Float64` arithmetic is used.

The precision is increased in two stages:

1. Newton's method evaluates the residual in double-double, but still solves with the `Float64` Jacobian and refines the solution iteratively.
2. If this refinement does not converge (it stagnates, diverges or gives non-finite values), the Jacobian, its LU factorization and the Taylor coefficients of the predictor are computed in double-double, too. This is called *full extended precision*. The tracker also enters it if the predictor needs double-double.

In full extended precision, each step additionally checks with a `Float64` step whether `Float64` suffices again. If so, the tracker leaves full extended precision. Otherwise the
path stays in extended precision. `init!` resets this state for each new path. The memory needed for full extended precision is allocated the first time it is used and reused for all later paths of the same tracker.

The Jacobian and the Taylor coefficients can only be computed in double-double for homotopies built from [`System`](@ref) or [`Homotopy`](@ref) (compiled, interpreted or
mixed), including parameter homotopies of such systems. For other [`AbstractHomotopy`](@ref)s, only the residual is evaluated in double-double.

Whether extended precision was used on a path can be read off the [`TrackerResult`](@ref), see the fields `extended_precision` and `extended_precision_used`.

### Paths that get stuck in extended precision

Close to a (numerically) singular point a path can get stuck: the steps are accepted, but they are so small, for example `1e-13`, that the path hardly advances, and neither extended precision nor smaller steps help. Such a path would only end after `max_steps` steps. 
To avoid this, the tracker terminates a path if, during `extended_precision_stall_steps` (default `100`) consecutive accepted steps in extended precision, the remaining distance to the target is reduced by less
than the fraction `extended_precision_stall_progress` (default `1e-8`). 
The tracking code is then `TrackerCode.terminated_step_size_too_small`. 
The distance is measured relative to the remaining distance, so
that paths approaching a singular endpoint geometrically are not affected. The check is disabled with `extended_precision_stall_steps = 0`.

