export is_at_infinity, is_failed

Base.@kwdef mutable struct EndgameTrackerOptions
    endgame_start::Float64 = 0.1
    terminate_ill_conditioned::Float64 = 1e12
    val_trust_tol::Float64 = 1e-3
    # singular solutions parameters
    val_singular_tol::Float64 = 1e-2
    min_cond_singular::Float64 = 1e6
    max_winding_number::Int = 20
    # at infinity parameters
    val_at_infinity_tol::Float64 = 1e-4
    strict_val_at_infinity_tol::Float64 = 1e-8
    min_cond_at_infinity::Float64 = 1e8
    max_t_at_infinity_without_cond::Float64 = 1e-14
    at_zero_check::Bool = false
    at_infinity_check::Bool = true
end

Base.show(io::IO, opts::EndgameTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::EndgameTrackerOptions) = opts

module EGTrackerReturnCode
import ..TrackerReturnCode
"""
    enum EGTrackerReturnCode

The possible states a `PathTracker` can be in:

* `EGTrackerReturnCode.tracking`
* `EGTrackerReturnCode.success`
* `EGTrackerReturnCode.at_infinity`
* `EGTrackerReturnCode.at_zero`
* `EGTrackerReturnCode.excess_solution`
* `EGTrackerReturnCode.post_check_failed`
* `EGTrackerReturnCode.terminated_accuracy_limit`
* `EGTrackerReturnCode.terminated_ill_conditioned`
* `EGTrackerReturnCode.terminated_invalid_startvalue`
* `EGTrackerReturnCode.terminated_max_winding_number`
* `EGTrackerReturnCode.terminated_max_iters`
* `EGTrackerReturnCode.terminated_step_size_too_small`
"""
@enum codes begin
    tracking
    success
    at_infinity
    at_zero
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_iters
    terminated_max_winding_number
    terminated_step_size_too_small
    terminated_unknown
    post_check_failed
    excess_solution
end

function Base.convert(::Type{codes}, code::TrackerReturnCode.codes)
    if code == TrackerReturnCode.success
        return success
    elseif code == TrackerReturnCode.terminated_max_iters
        return terminated_max_iters
    elseif code == TrackerReturnCode.terminated_accuracy_limit
        return terminated_accuracy_limit
    elseif code == TrackerReturnCode.terminated_ill_conditioned
        return terminated_ill_conditioned
    elseif code == TrackerReturnCode.terminated_invalid_startvalue
        return terminated_invalid_startvalue
    elseif code == TrackerReturnCode.terminated_step_size_too_small
        return terminated_step_size_too_small
    elseif code == TrackerReturnCode.terminated_unknown
        return terminated_unknown
    else
        return tracking
    end
end

end


"""
    is_success(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates a success in tracking.
"""
is_success(code::EGTrackerReturnCode.codes) = code == EGTrackerReturnCode.success

"""
    is_success(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates that a path diverged towards infinity.
"""
is_at_infinity(code::EGTrackerReturnCode.codes) = code == EGTrackerReturnCode.at_infinity

"""
    is_tracking(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates the tracking is not going on.
"""
is_tracking(code::EGTrackerReturnCode.codes) = code == EGTrackerReturnCode.tracking

"""
    is_invalid_startvalue(code::EGTrackerReturnCode.codes)

Returns `true` if the provided start value was not valid.
"""
is_invalid_startvalue(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.terminated_invalid_startvalue

"""
    is_terminated_callback(code::EGTrackerReturnCode.codes)

Returns `true` if the provided callback indicated a termination of the path.
"""
is_terminated_callback(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.terminated_callback


# State
Base.@kwdef mutable struct EndgameTrackerState{V<:AbstractVector}
    code::EGTrackerReturnCode.codes = EGTrackerReturnCode.tracking
    val::Valuation
    solution::V
    winding_number::Union{Nothing,Int} = nothing
    accuracy::Float64 = NaN
    prediction::V
    last_point::V
    last_t::Float64 = NaN
    max_winding_number_hit::Bool = false
    last_val_singular::Bool = false
    jump_to_zero_failed::Tuple{Bool,Bool} = (false, false)
end

EndgameTrackerState(x::AbstractVector) = EndgameTrackerState(;
    val = Valuation(length(x)),
    solution = copy(x),
    prediction = copy(x),
    last_point = copy(x),
)


struct EndgameTracker{
    H<:AbstractHomotopy,
    # V and V̄ need to have the same container type
    V<:AbstractVector{ComplexF64},
    V̄<:AbstractVector{ComplexDF64},
    M<:AbstractMatrix{ComplexF64},
}
    tracker::Tracker{H,V,V̄,M}
    state::EndgameTrackerState{V}
    options::EndgameTrackerOptions
end

function EndgameTracker(tracker::Tracker; kwargs...)
    options = EndgameTrackerOptions(; kwargs...)
    state = EndgameTrackerState(tracker.state.x)
    EndgameTracker(tracker, state, options)
end

Base.broadcastable(T::EndgameTracker) = Ref(T)

function init!(eg_tracker::EndgameTracker, x, t₁::Real)
    @unpack tracker, state, options = eg_tracker

    init!(tracker, x, t₁, 0.0)

    state.code = status(tracker)
    init!(state.val)
    state.solution .= NaN
    state.accuracy = NaN
    state.winding_number = nothing
    state.last_val_singular = false
    state.max_winding_number_hit = false
    state.jump_to_zero_failed = (false, false)

    eg_tracker
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
    cauchy!(p, tracker::CoreTracker, ::CauchyEndgame)

Try to predict the value of `x(0)` using the [`CauchyEndgame`](@ref).
For this we track the polygon defined by ``te^{i2πk/n}`` until we end again at ``x``.
Here ``n`` is the number of samples we take per loop, `samples_per_loop`.
The computation gives up if we have a winding number larger than `max_winding_number`.
It returns a tuple denoting the success ([`CauchyEndgameResult`](@ref)) the computed
winding number `m::Int` and th expected accuracy of the solution.

[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
function cauchy!(
    state::EndgameTrackerState,
    tracker::Tracker,
    options::EndgameTrackerOptions,
)
    @unpack last_point, prediction = state

    t = real(tracker.state.t)
    n₀ = max(ceil(Int, log(eps()) / log(t)), 3)
    @unpack x, μ, ω = tracker.state

    last_point .= x
    state.last_t = t
    prediction .= 0.0

    m = 1
    Δθ = 2π / n₀
    point_acc = tracker.state.accuracy
    result = CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    while m ≤ options.max_winding_number
        θⱼ = 0.0
        tⱼ = t
        for j = 1:n₀
            θⱼ += Δθ
            tⱼ = j == n₀ ? t : t * cis(θⱼ)
            res = track!(tracker, tⱼ)
            # TODO: Refine to guarantee high accuracy
            point_acc = max(point_acc, tracker.state.accuracy)

            if !is_success(res.code)
                result = CAUCHY_TERMINATED
                @goto _return
            end

            prediction .+= x
        end
        # Check that loop is closed
        if tracker.state.norm(last_point, x) < 4 * max(point_acc, tracker.state.accuracy)
            n = n₀ * m
            prediction .= prediction ./ n

            result = CAUCHY_SUCCESS
            break
        end

        m += 1
    end

    @label _return

    if result == CAUCHY_TERMINATED
        init!(tracker, last_point, t, 0.0, ω, μ; keep_steps = true)
    else
        init!(tracker, 0.0)
    end

    result, m, point_acc
end

function update_valuation!(state, tracker_state, t)
    @unpack x, x¹, x², x³, x⁴ = tracker_state
    update!(state.val, x, x¹, x², x³, x⁴, t)
end


function step!(eg_tracker::EndgameTracker, debug::Bool = false)
    @unpack tracker, state, options = eg_tracker

    proposed_t′ = real(tracker.state.t′)

    step!(tracker)

    state.code = status(tracker)

    if !is_tracking(state.code)
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x

        if state.code == EGTrackerReturnCode.terminated_accuracy_limit ||
           state.code == EGTrackerReturnCode.terminated_ill_conditioned

            verdict = analyze(
                state.val;
                finite_tol = options.val_trust_tol,
                singular_tol = options.val_singular_tol,
                max_winding_number = options.max_winding_number,
                at_infinity_tol = sqrt(options.val_at_infinity_tol),
                zero_is_finite = !options.at_zero_check,
            )

            if verdict.at_infinity
                return (state.code = EGTrackerReturnCode.at_infinity)
            end
        end
        return state.code
    end


    t = real(tracker.state.t)
    t < options.endgame_start || return state.code

    if tracker.state.last_step_failed
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), iszero(proposed_t′))
        return state.code
    else
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), false)
    end


    update_valuation!(state, tracker.state, t)

    if debug
        printstyled(
            "t = ",
            t,
            "  t′ = ",
            real(tracker.state.t′),
            "\n",
            color = tracker.state.extended_prec ? :blue : :yellow,
        )

        println("cond = ", LA.cond(tracker))
        println(state.val)
    end

    # analyze valuation
    verdict = analyze(
        state.val;
        finite_tol = options.val_trust_tol,
        singular_tol = options.val_singular_tol,
        max_winding_number = options.max_winding_number,
        at_infinity_tol = options.val_at_infinity_tol,
        zero_is_finite = !options.at_zero_check,
    )

    if debug
        println(verdict)
        println()
    end

    cond = LA.cond(tracker)

    verify_condition =
        cond > options.min_cond_at_infinity ||
        tracker.state.extended_prec || t < options.max_t_at_infinity_without_cond

    at_infinity =
        options.at_infinity_check &&
        (verdict.strict_at_infinity || (verdict.at_infinity && verify_condition))
    at_zero =
        options.at_zero_check &&
        (verdict.strict_at_zero || (verdict.at_zero && verify_condition))

    if at_infinity
        return (state.code = EGTrackerReturnCode.at_infinity)
    elseif at_zero
        return (state.code = EGTrackerReturnCode.at_zero)
    end

    # if valuation indicates singular solution
    # or finite but bad conditioned, start cauchy endgame
    use_cauchy_eg =
        (verdict.singular && verify_condition) || (
            verdict.finite && verdict.winding_number_candidate == 1 &&
            first(state.jump_to_zero_failed) && tracker.state.τ > 100t
        )

    # TODO: Check consistency of result -> second eg round
    if use_cauchy_eg # && state.last_val_singular
        res, m, acc_est = cauchy!(state, tracker, options)
        if res == CAUCHY_SUCCESS
            if m == verdict.winding_number_candidate
                state.winding_number = m
                state.solution .= state.prediction
                state.accuracy = acc_est
                state.code = EGTrackerReturnCode.success
            else
                state.max_winding_number_hit = true
                state.last_val_singular = false
            end
        elseif res == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
            if state.max_winding_number_hit
                state.code = EGTrackerReturnCode.terminated_max_winding_number
            else
                state.max_winding_number_hit = true
                state.last_val_singular = false
            end
        elseif res == CAUCHY_TERMINATED
            state.code = tracker.state.code
        end
    elseif use_cauchy_eg
        state.last_val_singular = true
    else
        state.last_val_singular = false
    end

    state.code
end

function track!(eg_tracker::EndgameTracker, x, t₁::Real; debug::Bool = false)
    init!(eg_tracker, x, t₁)

    while is_tracking(eg_tracker.state.code)
        step!(eg_tracker, debug)
    end

    eg_tracker.state.code
end

struct EGTrackerResult{V<:AbstractVector}
    return_code::EGTrackerReturnCode.codes
    solution::V
    t::Float64
    accuracy::Float64
    winding_number::Union{Nothing,Int}
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    valuation::Vector{Float64}
end

function EGTrackerResult(egtracker::EndgameTracker)
    @unpack tracker, state = egtracker
    t = real(tracker.state.t)
    EGTrackerResult(
        state.code,
        copy(state.solution),
        t,
        state.accuracy,
        state.winding_number,
        tracker.state.accepted_steps,
        tracker.state.rejected_steps,
        copy(state.val.val_x),
    )
end

Base.show(io::IO, r::EGTrackerResult) = print_fieldnames(io, r)
Base.show(io::IO, ::MIME"application/prs.juno.inline", r::EGTrackerResult) = r


function track(eg_tracker::EndgameTracker, x, t₁::Real = 1.0; debug::Bool = false)
    track!(eg_tracker, x, t₁; debug = debug)
    EGTrackerResult(eg_tracker)
end



"""
    solution(r::EGTrackerResult)

Get the solution of the path.
"""
solution(r::EGTrackerResult) = r.solution


"""
    accuracy(r::EGTrackerResult)

Get the accuracy of the solution. This is an estimate of the (relative) distance to the
true solution.
"""
accuracy(r::EGTrackerResult) = r.accuracy


"""
    winding_number(r::EGTrackerResult)

Get the winding number of the solution of the path. Returns `nothing` if it wasn't computed.
"""
winding_number(r::EGTrackerResult) = r.winding_number

"""
    is_success(r::EGTrackerResult)

Checks whether the path is successfull.
"""
is_success(r::EGTrackerResult) = is_success(r.return_code)

"""
    is_failed(r::EGTrackerResult)

Checks whether the path failed.
"""
is_failed(r::EGTrackerResult) = !(is_at_infinity(r) || is_success(r))

"""
    is_at_infinity(r::EGTrackerResult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::EGTrackerResult) = is_at_infinity(r.return_code)


"""
    is_finite(r::EGTrackerResult)

Checks whether the path result is finite.
"""
is_finite(r::EGTrackerResult) = is_success(r.return_code)
