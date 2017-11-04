## IETARATOR
function Base.start(tracker::Pathtracker)
    precondition!(tracker, tracker.low, tracker.low.cache)
    return 0
end

@inline function Base.next(tracker::Pathtracker, state)
    state += 1
    step!(tracker)

    tracker, state
end

@inline function Base.done(tracker::Pathtracker, state)
    if tracker.iter ≥ tracker.options.maxiters
        if tracker.options.verbose
            warn("Interrupted. Larger `maxiters` necessary.")
        end
        return true
    end
    if tracker.t ≤ 0.0
        return true
    end

    false
end
Base.eltype(::T) where {T<:Pathtracker} = T

#TODO: Should a user need more hooks?
@inline function step!(tracker::Pathtracker)
    @inbounds pre_perform_step!(tracker)

    # this is implemented from the different algorithms / caches
    if tracker.usehigh
        @inbounds perform_step!(tracker, tracker.high, tracker.high.cache)
    else
        @inbounds perform_step!(tracker, tracker.low, tracker.low.cache)
    end

    @inbounds post_perform_step!(tracker)
end

@inline function pre_perform_step!(tracker::Pathtracker)
    tracker.ds = tracker.steplength * tracker.sdiff
    tracker.snext = tracker.s + tracker.ds
    tracker.iter += 1
    nothing
end

@inline function post_perform_step!(tracker::Pathtracker)
    if tracker.step_sucessfull
        tracker.t -= tracker.steplength

        tracker.consecutive_successfull_steps += 1
        if tracker.consecutive_successfull_steps ==
           tracker.options.consecutive_successfull_steps_until_steplength_increase
            tracker.steplength *= tracker.options.steplength_increase_factor
            tracker.consecutive_successfull_steps = 0
        end

        tracker.s = tracker.snext

        copy_xnext_to_x!(tracker)
    else
        tracker.consecutive_successfull_steps = 0
        tracker.steplength *= tracker.options.steplength_decrease_factor
    end
    tracker.steplength = min(tracker.steplength, tracker.t)
    nothing
end

@inline function copy_xnext_to_x!(tracker::Pathtracker)
    if tracker.usehigh
        copy!(tracker.high.x, tracker.high.xnext)
    else
        copy!(tracker.low.x, tracker.low.xnext)
    end
    nothing
end

function run!(tracker::Pathtracker)
    precondition!(tracker, tracker.low, tracker.low.cache)
    while tracker.t > 0 && tracker.iter < tracker.options.maxiters
        step!(tracker)
    end
    if tracker.iter ≥ tracker.options.maxiters && tracker.options.verbose
        warn("Interrupted. Larger `maxiters` necessary.")
    end
end

function run!(
    tracker::Pathtracker{Low},
    startvalue::AbstractVector,
    start::Number=1.0,
    finish::Number=0.0) where {Low}
    reset!(tracker, startvalue, start, finish)
    run!(tracker)
    tracker
end
