module StepLength

export AbstractStepLength,
    AbstractStepLengthState,
    HeuristicStepLength,
    HeuristicStepLengthState,
    state,
    initial_steplength,
    update!

"""
    AbstractStepLength

Abstract type representing the configuration for the step length logic of a path tracker.
This always deals with a normalized parameter space, i.e., ``t ∈ [0, 1]``.
"""
abstract type AbstractStepLength end

"""
    AbstractStepLengthState

Abstract type representing the state of a `AbstractStepLength` implementation.
"""
abstract type AbstractStepLengthState end

"""
    state(::AbstractStepLength)::AbstractStepLengthState

Initialize the state for the given step length type.
"""
function state end

"""
    reset!(state::AbstractStepLengthState)

Reset the state to the initial state again.
"""
function reset! end


# FIXME: The API is not stable. This needs to be modified for verified pathtracking.

"""
    isrelative(::AbstractStepLength)

Indicate whether the returned step length is relative or absolute, i.e.,
if we want to track from `1` to `5`. A relative step length of `0.1` would result
in an effective step length of `0.4` while an absolute step length of `0.1` would result
in an effective step length of `0.1`.
"""
function isrelative end

"""
    initial_steplength(::AbstractStepLength)

Return the initial step length of a path.
"""
function initial_steplength end

"""
    update(curr_steplength, step::AbstractStepLength, state::AbstractStepLengthState, success::Bool)

Returns `(t, status)` where `t` is the new step length and `status` is either
`:ok` (usually) or  any symbol indicating why the update failed (e.g. `:steplength_too_small` if the step size is too small).
The argument `success` indicates whether the step was successfull. This
modifies `state` if necessary.
"""
function update end


# Implementation
"""
    HeuristicStepLength(;initial=0.1,
        increase_factor=2.0,
        decrease_factor=inv(increase_factor),
        consecutive_successes_necessary=3,
        maximal_steplength=0.1
        minimal_steplength=1e-14)

The step length is defined as follows. Initially the step length is `initial`.
If `consecutive_successes_necessary` consecutive steps were sucessfull the step length
is increased by the factor `increase_factor`. If a step fails, i.e. the corrector does
not converge, the steplength is reduced by the factor `decrease_factor`.
"""
struct HeuristicStepLength <: AbstractStepLength
    initial::Float64
    increase_factor::Float64
    decrease_factor::Float64
    consecutive_successes_necessary::Int
    maximal_steplength::Float64
    minimal_steplength::Float64

end
function HeuristicStepLength(;initial=0.1,
    increase_factor=2.0,
    decrease_factor=inv(increase_factor),
    consecutive_successes_necessary=5,
    maximal_steplength=0.1,
    minimal_steplength=1e-14)

    HeuristicStepLength(initial, increase_factor,
        decrease_factor, consecutive_successes_necessary, maximal_steplength,
        minimal_steplength)
end

mutable struct HeuristicStepLengthState <: AbstractStepLengthState
    consecutive_successes::Int
end

state(step::HeuristicStepLength) = HeuristicStepLengthState(0)
function reset!(state::HeuristicStepLengthState)
    state.consecutive_successes = 0
end

isrelative(::HeuristicStepLength) = true
initial_steplength(step::HeuristicStepLength) = step.initial

function update(curr_steplength, step::HeuristicStepLength, state::HeuristicStepLengthState, success)
    if success
        if (state.consecutive_successes += 1) == step.consecutive_successes_necessary
            state.consecutive_successes = 0
            new_steplength = min(step.maximal_steplength, step.increase_factor * curr_steplength)
            return (new_steplength, :ok)
        end
        return (curr_steplength, :ok)
    end
    # reset successes
    state.consecutive_successes = 0

    # we decrease the steplength
    new_steplength = step.decrease_factor * curr_steplength

    if new_steplength < step.minimal_steplength
        return (new_steplength, :steplength_too_small)
    end

    (new_steplength, :ok)
end

end
