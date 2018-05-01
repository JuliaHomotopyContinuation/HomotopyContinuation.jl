module StepLength

export AbstractStepLength,
    AbstractStepLengthState,
    HeuristicStepLength,
    HeuristicStepLengthState,
    state,
    reset!,
    relsteplength,
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
    state(::AbstractStepLength, start, target)::AbstractStepLengthState

Initialize the state for the given step length type.
"""
function state end

"""
    reset!(state::AbstractStepLengthState, step::AbstractStepLength, start, target)

Reset the state to the initial state again.
"""
function reset! end


# FIXME: The API is not stable. This needs to be modified for verified pathtracking.
"""
    current_steplength(::AbstractStepLengthState)

Return the current step length of a path.
"""
function current_steplength end

"""
    relsteplength(state)

Return the current steplength in a relative measure. This is always a positive
number and always between 0.0 and 1.0.
"""
function relsteplength end

"""
    update!(state::AbstractStepLengthState, step::AbstractStepLength, success::Bool)

Returns a status which is either
`:ok` (usually) or  any symbol indicating why the update failed (e.g. `:steplength_too_small` if the step size is too small).
The argument `success` indicates whether the step was successfull.
"""
function update! end


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
    maximal_steplength=max(0.1, initial),
    minimal_steplength=1e-14)

    HeuristicStepLength(initial, increase_factor,
        decrease_factor, consecutive_successes_necessary, maximal_steplength,
        minimal_steplength)
end

mutable struct HeuristicStepLengthState <: AbstractStepLengthState
    consecutive_successes::Int
    steplength::Float64
    length_start_target::Float64
end

function state(step::HeuristicStepLength, start, target)
    length_start_target = convert(Float64, abs(start-target))
    steplength = min(step.initial, length_start_target)
    HeuristicStepLengthState(0, steplength, length_start_target)
end
function reset!(state::HeuristicStepLengthState, step::HeuristicStepLength, start, target)
    state.consecutive_successes = 0
    state.length_start_target = convert(Float64, abs(start-target))
    state.steplength = min(step.initial, state.length_start_target)
end


relsteplength(state::HeuristicStepLengthState) = min(state.steplength / state.length_start_target, 1.0)

function update!(state::HeuristicStepLengthState, step::HeuristicStepLength, success)
    curr_steplength = state.steplength
    if success
        if (state.consecutive_successes += 1) == step.consecutive_successes_necessary
            # try to increase step length
            state.consecutive_successes = 0
            state.steplength = min(
                step.maximal_steplength,
                step.increase_factor * curr_steplength,
                state.length_start_target)
            return :ok
        end
        return :ok
    end
    # reset successes
    state.consecutive_successes = 0

    # we decrease the steplength
    state.steplength = step.decrease_factor * curr_steplength

    if state.steplength < step.minimal_steplength
        return :steplength_too_small
    end

    :ok
end

end
