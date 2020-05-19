export PathTracker, PathTrackerOptions, track

"""
    AbstractPathTracker

Supertype for path trackers.
"""
abstract type AbstractPathTracker end

Base.broadcastable(T::AbstractPathTracker) = Ref(T)

"""
    PathTrackerOptions(; options...)

Options controlling the behaviour of a [`PathTracker`](@ref).

## Options

* `endgame_start::Float64 = 0.1`: The point `t` in time where the endgame starts.

### Endgame parameters
These parameters control the behaviour during the endgame. See [^BT20] for details.

* `min_cond::Float64 = 1e6`: The minimal growth of the condition number after which an
  endgame strategy is considered to be applied.
* `min_growth::Float64 = 100`: The minimal relative growth of a coordinate necessary to
  to be considered going to infininity (resp. zero).
* `at_infinity_check::Bool = true`: Whether divering paths should be truncated.
* `zero_is_at_infinity::Bool = false`: Whether paths going to a solution where at least one
  coordinates is zero should also be considered diverging.
* `val_finite_tol::Float64 = 1e-2`: Tolerance on the valuation which has to be satisfied
  before the Cauchy endgame is started.
* `max_winding_number::Int = 12`: The maximal winding number which is attempted in the
  Cauchy endgame.
* `val_at_infinity_tol::Float64 = 1e-3`: Tolerance on the valuation which has to be
  satisfied before a path is considered to diverge / go to infinity.

[^BT20]: Breiding, P. and Timme, S. "Tropical Endgame", In preparation (2020)
"""
Base.@kwdef mutable struct PathTrackerOptions
    endgame_start::Float64 = 0.1
    max_endgame_steps::Int = 2000
    # eg parameter
    min_cond::Float64 = 1e4
    min_coord_growth::Float64 = 100.0
    zero_is_at_infinity::Bool = false
    at_infinity_check::Bool = true
    # valuation etc
    val_finite_tol::Float64 = 1e-3
    val_at_infinity_tol::Float64 = 1e-3

    # singular solutions parameters
    max_winding_number::Int = 12
end

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts

"""
    PathTrackerCode

The possible states a `PathTracker` can be in:

* `PathTrackerCode.tracking`
* `PathTrackerCode.success`
* `PathTrackerCode.at_infinity`
* `PathTrackerCode.at_zero`
* `PathTrackerCode.excess_solution`
* `PathTrackerCode.post_check_failed`
* `PathTrackerCode.polyhedral_failed`
* `PathTrackerCode.terminated_accuracy_limit`
* `PathTrackerCode.terminated_ill_conditioned`
* `PathTrackerCode.terminated_invalid_startvalue`
* `PathTrackerCode.terminated_max_winding_number`
* `PathTrackerCode.terminated_max_steps`
* `PathTrackerCode.terminated_step_size_too_small`
"""
module PathTrackerCode
import ..TrackerCode

@enum codes begin
    tracking
    success
    at_infinity
    at_zero
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_steps
    terminated_max_winding_number
    terminated_step_size_too_small
    terminated_unknown
    post_check_failed
    excess_solution
    polyhedral_failed
end

function Base.convert(::Type{codes}, code::TrackerCode.codes)
    if code == TrackerCode.success
        return success
    elseif code == TrackerCode.terminated_max_steps
        return terminated_max_steps
    elseif code == TrackerCode.terminated_accuracy_limit
        return terminated_accuracy_limit
    elseif code == TrackerCode.terminated_ill_conditioned
        return terminated_ill_conditioned
    elseif code == TrackerCode.terminated_invalid_startvalue
        return terminated_invalid_startvalue
    elseif code == TrackerCode.terminated_step_size_too_small
        return terminated_step_size_too_small
    elseif code == TrackerCode.terminated_unknown
        return terminated_unknown
    else
        return tracking
    end
end

end


"""
    is_success(code::PathTrackerCode.codes)

Returns `true` if `status` indicates a success in tracking.
"""
is_success(code::PathTrackerCode.codes) = code == PathTrackerCode.success

"""
    is_success(code::PathTrackerCode.codes)

Returns `true` if `status` indicates that a path diverged towards infinity.
"""
is_at_infinity(code::PathTrackerCode.codes) = code == PathTrackerCode.at_infinity

"""
    is_tracking(code::PathTrackerCode.codes)

Returns `true` if `status` indicates the tracking is not going on.
"""
is_tracking(code::PathTrackerCode.codes) = code == PathTrackerCode.tracking

"""
    is_invalid_startvalue(code::PathTrackerCode.codes)

Returns `true` if the provided start value was not valid.
"""
is_invalid_startvalue(code::PathTrackerCode.codes) =
    code == PathTrackerCode.terminated_invalid_startvalue

###
### State
###
Base.@kwdef mutable struct PathTrackerState
    code::PathTrackerCode.codes = PathTrackerCode.tracking
    val::Valuation
    endgame_started::Bool = false
    row_scaling::Vector{Float64}
    col_scaling::Vector{Float64}
    cond_eg_start::Float64 = 1.0
    steps_eg::Int = 0
    solution::Vector{ComplexF64}
    t_last_cauchy::Float64 = NaN
    winding_number::Union{Nothing,Int} = nothing
    accuracy::Float64 = NaN
    cond::Float64 = 1.0
    prediction::Vector{ComplexF64}
    last_point::Vector{ComplexF64}
    last_t::Float64 = NaN
    max_winding_number_hit::Bool = false
    cauchy_failures::Int = 0
    jump_to_zero_failed::Tuple{Bool,Bool} = (false, false)
end

PathTrackerState(npolynomials::Integer, x::AbstractVector) = PathTrackerState(;
    val = Valuation(length(x)),
    row_scaling = zeros(npolynomials),
    col_scaling = zeros(length(x)),
    solution = copy(x),
    prediction = copy(x),
    last_point = copy(x),
)

"""
    PathTracker(tracker::Tracker; options = PathTrackerOptions())
    PathTracker(H::AbstractHomotopy; options = PathTrackerOptions())

A `PathTracker` combines a [`Tracker`](@ref) with an endgame. That is,
while a [`Tracker`](@ref) assumes that the solution path is non-singular and convergent, the endgame
allows to handle singular endpoints as well as diverging paths.
To compute singular solutions the *Cauchy endgame* used, for divering paths a strategy
based on the valuation of local Puiseux series expansion of the path is used.
See [^BT20] for a detailed description.
By convention, a `PathTracker` always tracks from ``t=1`` to ``t = 0``.
See [`PathTrackerOptions`](@ref) for the possible options.

[^BT20]: Breiding, P. and Timme, S. "Tropical Endgame", In preparation (2020)
"""
struct PathTracker{
    H<:AbstractHomotopy,
    N, # AutomaticDifferentiation
    M<:AbstractMatrix{ComplexF64},
} <: AbstractPathTracker
    tracker::Tracker{H,N,M}
    state::PathTrackerState
    options::PathTrackerOptions
end

PathTracker(H::Homotopy; kwargs...) = PathTracker(ModelKitHomotopy(H); kwargs...)
function PathTracker(H::AbstractHomotopy; tracker_options = TrackerOptions(), kwargs...)
    PathTracker(Tracker(H; options = tracker_options); kwargs...)
end
function PathTracker(tracker::Tracker; options = PathTrackerOptions())
    state = PathTrackerState(size(tracker.homotopy, 1), tracker.state.x)
    PathTracker(tracker, state, options)
end

Base.broadcastable(T::PathTracker) = Ref(T)

function init!(path_tracker::PathTracker, x, t₁::Real; ω::Float64 = NaN, μ::Float64 = NaN)
    @unpack tracker, state, options = path_tracker

    init!(tracker, x, t₁, 0.0; ω = ω, μ = μ)

    state.code = status(tracker)
    init!(state.val)
    state.endgame_started = false
    state.row_scaling .= 1
    state.col_scaling .= weights(tracker.state.norm)
    state.cond_eg_start = 1.0
    state.cond = 1.0
    state.solution .= NaN
    state.accuracy = NaN
    state.t_last_cauchy = NaN
    state.winding_number = nothing
    state.steps_eg = 0
    state.max_winding_number_hit = false
    state.cauchy_failures = 0
    state.jump_to_zero_failed = (false, false)

    path_tracker
end

"""
    CauchyEndgameResult

An enum indicating the result of the [`cauchy!`](@ref) computation.

# Cases
* `CAUCHY_SUCCESS`: The endgame was successfull.
* `CAUCHY_TERMINATED_MAX_WINDING_NUMBER`: The endgame was terminated since the winding
  number is larger than the provided threshold.
* `CAUCHY_TERMINATED`: The endgame was terminated due to some other error in the path
  tracking.
"""
@enum CauchyEndgameResult begin
    CAUCHY_SUCCESS
    CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    CAUCHY_TERMINATED
end

"""
    cauchy!(state::PathTrackerState, tracker::Tracker, options::PathTrackerOptions)

Try to predict the value of `x(0)` using the [`CauchyEndgame`](@ref).
For this we track the polygon defined by ``te^{i2πk/n}`` until we end again at ``x``.
Here ``n`` is the number of samples we take per loop, `samples_per_loop`.
The computation gives up if we have a winding number larger than `max_winding_number`.
It returns a tuple denoting the success ([`CauchyEndgameResult`](@ref)) the computed
winding number `m::Int` and th expected accuracy of the solution.

[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
function cauchy!(state::PathTrackerState, tracker::Tracker, options::PathTrackerOptions)
    @unpack last_point, prediction = state

    t = real(tracker.state.t)
    # Mathemtically, we only need `n₀ = ceil(Int, log(eps()) / log(t))` many sample points
    # to achieve the maximal accuracy since the error is ≈ t^n₀.
    # However, we have to be careful that we are not getting too close to the singularity
    # during tracking. E.g. with `n₀ = 3` we track fairly close to the origin for the
    # first first and third path. So we require at least 8 sample points.
    n₀ = max(ceil(Int, log(eps()) / log(t)), 8)
    @unpack x, μ, ω = tracker.state

    # always use extended precision for cauchy endgame
    prediction_acc = use_extended_precision!(tracker)
    # fix tracker to not flip between extended precision and and mach. precision
    tracker.state.keep_extended_prec = true
    state.last_point .= tracker.state.x
    state.last_t = tracker.state.t
    prediction .= 0.0
    sample_point_acc = Inf
    m = 1
    Δθ = 2π / n₀
    result = CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    while m ≤ options.max_winding_number
        θⱼ = 0.0
        tⱼ = t
        for j = 1:n₀
            θⱼ += Δθ
            tⱼ = j == n₀ ? t : t * cis(θⱼ)
            res = track!(tracker, tⱼ)
            sample_point_acc = tracker.state.accuracy
            prediction_acc = max(prediction_acc, sample_point_acc)

            if !is_success(res)
                result = CAUCHY_TERMINATED
                @goto _return
            end

            prediction .+= x
        end
        # Check that loop is closed
        d = tracker.state.norm(last_point, x)
        if d < 10 * max(prediction_acc, sample_point_acc)
            n = n₀ * m
            prediction .= prediction ./ n

            result = CAUCHY_SUCCESS
            break
        end
        m += 1
    end

    @label _return

    if result == CAUCHY_TERMINATED
        init!(tracker, last_point, t, 0.0; ω = ω, μ = μ, keep_steps = true)
    else
        init!(tracker, 0.0)
    end

    result, m, prediction_acc
end

function max_ratio(u::Vector{Float64}, v::Vector{Float64})
    max_rat = -Inf
    for i in eachindex(u)
        max_rat = max(max_rat, u[i] / v[i])
    end
    max_rat
end
function min_ratio(u::Vector{Float64}, v::Vector{Float64})
    min_rat = Inf
    for i in eachindex(u)
        min_rat = max(min_rat, u[i] / v[i])
    end
    min_rat
end

function step!(path_tracker::PathTracker, debug::Bool = false)
    @unpack tracker, state, options = path_tracker

    proposed_t′ = real(tracker.state.t′)
    state.last_point .= tracker.state.x
    state.last_t = tracker.state.t
    step!(tracker, debug)

    state.code = status(tracker)

    if !is_tracking(state.code)
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        state.winding_number = nothing
        state.cond = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)

        # If a path got terminated, let's be more generous to classify them as at_infinity
        if state.code == PathTrackerCode.terminated_accuracy_limit ||
           state.code == PathTrackerCode.terminated_ill_conditioned

            @label terminated
            res = analyze(
                state.val;
                finite_tol = 0.0,
                # at_infinity_tol = sqrt(options.val_at_infinity_tol),
                zero_is_finite = !options.zero_is_at_infinity,
            )


            if options.at_infinity_check &&
               res.at_infinity_tol < cbrt(options.val_at_infinity_tol)
                return (state.code = PathTrackerCode.at_infinity)
            elseif options.zero_is_at_infinity &&
                   res.at_zero_tol < cbrt(options.val_at_infinity_tol)
                return (state.code = PathTrackerCode.at_zero)
            end
        end
        return state.code
    end


    t = real(tracker.state.t)
    t < options.endgame_start || return state.code

    if !state.endgame_started
        #=
         We enter the endgame.
         We monitor the condition number with a fixed scaling of the columns and rows,
         which is recorded now.
        =#
        state.col_scaling .= weights(tracker.state.norm)
        row_scaling!(state.row_scaling, jacobian(tracker.state.jacobian), state.col_scaling)
        κ = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
        state.cond_eg_start = κ

        tracker.state.use_strict_β_τ = true
        state.endgame_started = true
    end

    state.steps_eg += 1

    if state.steps_eg > options.max_endgame_steps
        return (state.code = PathTrackerCode.terminated_max_steps)
    end
    if tracker.state.last_step_failed
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), iszero(proposed_t′))
        return state.code
    else
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), false)
    end


    update!(state.val, tracker.predictor, t)
    (val_finite, at_infinity_tol, at_zero_tol) = analyze(
        state.val;
        finite_tol = options.val_finite_tol * exp10(-state.cauchy_failures),
        # at_infinity_tol = options.val_at_infinity_tol,
        zero_is_finite = !options.zero_is_at_infinity,
    )


    if val_finite || (min(at_zero_tol, at_infinity_tol) < sqrt(options.val_at_infinity_tol))
        state.cond = egcond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
    end

    # For the truncation of paths use the stronger requirement of a relative
    # condition number increase
    κ = state.cond / state.cond_eg_start

    κ_min = options.min_cond

    val_at_infinity_tol = begin
        if κ > κ_min^3
            cbrt(options.val_at_infinity_tol)
        elseif κ > κ_min^2
            sqrt(options.val_at_infinity_tol)
        elseif κ > κ_min
            options.val_at_infinity_tol
        else
            0.0
        end
    end

    at_infinity =
        options.at_infinity_check &&
        at_infinity_tol < val_at_infinity_tol && (
            validate_coord_growth(
                state.val,
                state.col_scaling,
                weights(tracker.state.norm);
                finite_tol = options.val_finite_tol,
                at_infinity_tol = val_at_infinity_tol,
                tol = options.min_coord_growth,
            ) || tracker.state.extended_prec
        )
    at_zero =
        options.at_infinity_check &&
        options.zero_is_at_infinity &&
        at_infinity_tol < val_at_infinity_tol && (
            validate_coord_growth(
                state.val,
                weights(tracker.state.norm),
                state.col_scaling;
                finite_tol = options.val_finite_tol,
                at_infinity_tol = val_at_infinity_tol,
                tol = options.min_coord_growth,
            ) || tracker.state.extended_prec
        )
    finite = val_finite && κ > options.min_cond

    if debug
        color = tracker.state.extended_prec ? :blue : :yellow
        printstyled("t = ", t, "  t′ = ", real(tracker.state.t′), "\n", color = color)
        κ = egcond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
        @show tracker.predictor.cond_H_ẋ
        print("κ = ")
        Printf.@printf("%.3e (%.3e)\n", κ, κ / state.cond_eg_start)
        println(state.val)

        printstyled("Val: ", bold = true)
        if val_finite
            printstyled("val_finite ", color = :blue, bold = true)
        end
        if finite
            printstyled("finite ", color = :blue, bold = true)
        end
        if at_infinity
            printstyled("at_infinity ", color = :blue, bold = true)
        end
        if at_zero
            printstyled("at_zero ", color = :blue, bold = true)
        end
        println("\n")
    end

    if at_infinity
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        return (state.code = PathTrackerCode.at_infinity)

    elseif at_zero
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        return (state.code = PathTrackerCode.at_zero)

    elseif finite && (isnan(state.t_last_cauchy) || 10t < state.t_last_cauchy)
        res, m, acc_est = cauchy!(state, tracker, options)
        if debug
            printstyled("Cauchy result: ", res, " ", m, " ", acc_est, "\n"; color = :blue)
        end
        @show state.prediction
        if res == CAUCHY_SUCCESS
            if state.winding_number === nothing
                @label save_cauchy_result
                state.t_last_cauchy = t
                state.winding_number = m
                state.solution .= state.prediction
                state.accuracy = acc_est
            elseif state.winding_number != m
                state.cauchy_failures += 1
                @goto save_cauchy_result
            elseif state.winding_number == m
                d =
                    distance(state.prediction, state.solution, InfNorm()) /
                    maximum(weights(tracker.state.norm))
                @show d, state.accuracy, acc_est
                if d < 100 * max(state.accuracy, acc_est)
                    state.solution .= state.prediction
                    state.accuracy = d
                    state.cond = LA.cond(
                        tracker,
                        state.solution,
                        0.0,
                        state.row_scaling,
                        state.col_scaling,
                    )
                    return (state.code = PathTrackerCode.success)
                else
                    state.accuracy = acc_est
                    state.solution .= state.prediction
                    state.cauchy_failures += 1
                    @goto save_cauchy_result
                end
            end
        elseif res == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
            if state.max_winding_number_hit
                state.code = PathTrackerCode.terminated_max_winding_number
            else
                state.cauchy_failures += 1
                state.max_winding_number_hit = true
            end
        elseif res == CAUCHY_TERMINATED
            state.code = tracker.state.code
        end
    end

    # Catch ill behavior and terminate the tracking
    # 1) check cond(H_x, ẋ)
    if tracker.predictor.cond_H_ẋ > 1e13 && maximum(state.val.val_tẋ) > 1e2
        state.code = PathTrackerCode.terminated_ill_conditioned
        @goto terminated
    end

    state.code
end


function track!(
    path_tracker::PathTracker,
    x,
    t₁::Real;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    debug::Bool = false,
)
    init!(path_tracker, x, t₁; ω = ω, μ = μ)

    while is_tracking(path_tracker.state.code)
        step!(path_tracker, debug)
    end

    path_tracker.state.code
end

function solution(path_tracker::PathTracker)
    get_solution(path_tracker.tracker.homotopy, path_tracker.state.solution, 0.0)
end


function PathResult(
    path_tracker::PathTracker,
    start_solution = nothing,
    path_number = nothing,
)
    @unpack tracker, state, options = path_tracker
    H = tracker.homotopy


    if is_success(state.code)
        t = 0.0
        solution = get_solution(H, state.solution, 0.0)
        evaluate!(tracker.corrector.r, H, state.solution, complex(0.0))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    else
        t = real(tracker.state.t)
        solution = get_solution(H, tracker.state.x, t)
        evaluate!(tracker.corrector.r, H, tracker.state.x, complex(t))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    end

    PathResult(
        return_code = Symbol(state.code),
        solution = solution,
        t = t,
        accuracy = state.accuracy,
        residual = residual,
        condition_jacobian = state.cond,
        winding_number = state.winding_number,
        multiplicity = nothing,
        last_path_point = (get_solution(H, state.last_point, state.last_t), state.last_t),
        valuation = t > options.endgame_start ? nothing : copy(state.val.val_x),
        start_solution = start_solution,
        path_number = path_number,
        ω = tracker.state.ω,
        μ = tracker.state.μ,
        extended_precision = tracker.state.extended_prec,
        accepted_steps = tracker.state.accepted_steps,
        rejected_steps = tracker.state.rejected_steps,
        steps_eg = state.steps_eg,
        extended_precision_used = tracker.state.used_extended_prec,
    )
end

"""
    track(path_tracker::PathTracker, x::AbstractVector, t::Real = 1.0;
          path_number = nothing, debug = false)

Track the given start solution `x` from `t` towards `0` using the given `path_tracker`.
Returns a [`PathResult`](@ref).

    track(path_tracker::PathTracker, r::PathResult, t::Real = 1.0;
          path_number = nothing, debug = false)

Track `solution(r)` from `t` towards `0` using the given `path_tracker`.
"""
function track(
    path_tracker::PathTracker,
    x::AbstractVector,
    t₁::Real = 1.0;
    path_number::Union{Nothing,Int} = nothing,
    kwargs...,
)
    track!(path_tracker, x, t₁; kwargs...)
    PathResult(path_tracker, x, path_number)
end
function track(path_tracker::PathTracker, r::PathResult, t₁::Real; kwargs...)
    track(path_tracker, solution(r), t₁; ω = r.ω, μ = r.μ, kwargs...)
end
