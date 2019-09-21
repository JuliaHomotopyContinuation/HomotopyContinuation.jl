"""
    enum PathTrackerStatus

The possible states a [`PathTracker`](@ref) can be in:
* `PT_TRACKING`
* `PT_SUCCESS`
* `PT_AT_INFINITY`
* `PT_EXCESS_SOLUTION`
* `PT_POST_CHECK_FAILED`
* `PT_TERMINATED_ACCURACY_LIMIT`
* `PT_TERMINATED_ILL_CONDITIONED`
* `PT_TERMINATED_INVALID_STARTVALUE`
* `PT_TERMINATED_MAX_ITERS`
* `PT_TERMINATED_STEP_SIZE_TOO_SMALL`
"""
@enum PathTrackerStatus begin
    PT_TRACKING
    PT_SUCCESS
    PT_AT_INFINITY
    PT_TERMINATED_MAX_ITERS
    PT_TERMINATED_INVALID_STARTVALUE
    PT_TERMINATED_STEP_SIZE_TOO_SMALL
    PT_TERMINATED_ACCURACY_LIMIT
    PT_TERMINATED_ILL_CONDITIONED
    PT_POST_CHECK_FAILED
    PT_EXCESS_SOLUTION
end

"""
    path_tracker_status(code::CoreTrackerStatus)
Construct a [`PathTrackerStatus`](@ref) from a [`CoreTrackerStatus`](@ref).
"""
function path_tracker_status(code::CoreTrackerStatus)
    if code == CT_SUCCESS
        return PT_SUCCESS
    elseif code == CT_TERMINATED_INVALID_STARTVALUE
        return PT_TERMINATED_INVALID_STARTVALUE
    elseif code == CT_TERMINATED_MAX_ITERS
        return PT_TERMINATED_MAX_ITERS
    elseif code == CT_TERMINATED_STEP_SIZE_TOO_SMALL
        return PT_TERMINATED_STEP_SIZE_TOO_SMALL
    elseif code == CT_TERMINATED_ILL_CONDITIONED
        return PT_TERMINATED_ILL_CONDITIONED
    else
        return PT_TRACKING
    end
end

is_success(status::PathTrackerStatus) = status == PT_SUCCESS
is_tracking(status::PathTrackerStatus) = status == PT_TRACKING

mutable struct PathTrackerState{AV<:AbstractVector}
    status::PathTrackerStatus
    valuation::Valuation
    prediction::AV
    winding_number::Union{Nothing,Int}
    s::Float64
end

function PathTrackerState(x::AbstractVector; at_infinity_check::Bool = true)
    status = PT_TRACKING
    valuation = Valuation(x; affine = at_infinity_check)
    prediction = copy(x)
    winding_number = nothing
    s = 0.0
    PathTrackerState(status, valuation, prediction, winding_number, s)
end

function init!(state::PathTrackerState, s::Float64)
    state.status = PT_TRACKING
    init!(state.valuation)
    state.prediction .= 0.0
    state.winding_number = nothing
    state.s = s
end

mutable struct PathTrackerOptions
    endgame_start::Float64
    at_infinity_check::Bool
end

struct PathTracker{
    AV<:AbstractVector{Complex{Float64}},
    Prob<:AbstractProblem,
    CT<:CoreTracker{Float64,AV},
}
    problem::Prob
    core_tracker::CT
    endgame::CauchyEndgame{AV}
    state::PathTrackerState{AV}
    options::PathTrackerOptions
end



function PathTracker(
    prob::Problem,
    core_tracker::CoreTracker;
    at_infinity_check = default_at_infinity_check(prob),
    samples_per_loop::Int = 8,
    max_winding_number::Int = 12,
    endgame_start::Float64 = 2.0,
)
    state = PathTrackerState(current_x(core_tracker); at_infinity_check = at_infinity_check)
    endgame = CauchyEndgame(
        current_x(core_tracker);
        samples_per_loop = samples_per_loop,
        max_winding_number = max_winding_number,
    )
    options = PathTrackerOptions(endgame_start, at_infinity_check)

    PathTracker(prob, core_tracker, endgame, state, options)
end

default_at_infinity_check(prob::Problem{AffineTracking}) = true
default_at_infinity_check(prob::Problem{ProjectiveTracking}) = homvars(prob) !== nothing

function track!(
    tracker::PathTracker,
    x,
    s₁ = 0.0,
    s₀ = -log(10 * tracker.core_tracker.options.min_step_size),
)
    init!(tracker, x, s₁, s₀)
    while is_tracking(tracker.state.status)
        step!(tracker)
    end
    tracker.state.status
end

function init!(
    tracker::PathTracker,
    x,
    s₁ = 0.0,
    s₀ = -log(10 * tracker.core_tracker.options.min_step_size),
)
    init!(tracker.core_tracker, x, s₁, s₀)
    init!(tracker.state, s₁)
    tracker
end

function step!(tracker::PathTracker)
    @unpack core_tracker, state, options, endgame = tracker

    step!(core_tracker)
    # If we didn't move forward there is nothing to do
    core_tracker.state.last_step_failed && return state.status
    state.s = real(current_t(core_tracker))

    ct_status = status(core_tracker)
    if is_tracking(ct_status)
        # If we are not yet in the region where we consider endgames move on
        state.s < options.endgame_start && return state.status

        # update the valuation
        update!(state.valuation, core_tracker)

        # We only care about the valuation judgement if the path is starting to get worse
        sign_of_ill_conditioning(core_tracker; safety_factor = 1e4) || return state.status

        # Judge the current valuation to determine how to proceed
        verdict = judge(state.valuation; tol = 1e-3, tol_at_infinity = 1e-4)
        if options.at_infinity_check && verdict == VAL_AT_INFINITY
            status.status = PT_AT_INFINITY
        # If we expect a finite value let's do the cauchy endgame
        elseif verdict == VAL_FINITE
            # Perform endgame to estimate singular solution
            retcode, m = predict!(state.prediction, core_tracker, endgame)
            if retcode == CAUCHY_SUCCESS
                state.winding_number = m
                state.status = PT_SUCCESS
            elseif retcode == CAUCHY_TERMINATED_ACCURACY_LIMIT
                # TODO: Allow adaptive precision depending on the accuracy?
                state.status = PT_TERMINATED_ILL_CONDITIONED
            elseif retcode == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
                # TODO
            else
                state.status = PT_TERMINATED_ILL_CONDITIONED
            end
        end
    elseif is_success(ct_status)
        # TODO: We shouldn't declare immediately success, check that we are actually there
        # otherwise we need to move further (and maybe decrease min_step_size)
        state.status = PT_SUCCESS

    elseif ct_status == CT_TERMINATED_ACCURACY_LIMIT
        # TODO: What exactly do here?
        state.status = PT_TERMINATED_ACCURACY_LIMIT
    end

    state.status
end

"Update the valuation with the current state."
update!(val::Valuation, T::CoreTracker) =
    update!(val, T.state.x, T.state.ẋ, real(current_t(T)), T.predictor)

"""
    sign_of_ill_conditioning(CT::CoreTracker; safety_factor::Float64 = 1e4)

Returns `true` if the path has some signs of ill-conditioning. For this we check whether
the condition number of the Jacboan is larger than `safety_factor` or if
the safety factor to reaching limit accuracy is less than `safety_factor`.
"""
function sign_of_ill_conditioning(CT::CoreTracker; safety_factor::Float64 = 1e4)
    cond(CT) > safety_factor || safety_factor * CT.state.limit_accuracy > CT.options.accuracy
end

function Base.iterate(tracker::PathTracker, state::Int = 0)
    if is_tracking(tracker.state.status)
        step!(tracker)
        tracker, state + 1
    else
        nothing
    end
end
