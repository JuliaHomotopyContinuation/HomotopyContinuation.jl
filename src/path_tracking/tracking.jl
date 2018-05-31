"""
    track(tracker, x₁, t₁, t₀)::PathTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies tracker.

    track(tracker, xs, t₁, t₀)::Vector{PathTrackerResult}

Track all values in xs from `t₁` to `t₀` using the given `PathTracker` `tracker`.
"""
function track(tracker::PathTracker, x₁::AbstractVector, t₁, t₀; kwargs...)
     track!(tracker, x₁, t₁, t₀; kwargs...)
     PathTrackerResult(tracker)
end

"""
     track!(x₀, tracker, x₁, t₁, t₀)

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`
and store the result in `x₀` if the tracking was successfull.
Returns a `Symbol` indicating the status.
If the tracking was successfull it is `:success`.
"""
function track!(x₀, tracker::PathTracker, x₁, t₁, t₀)
     track!(tracker, x₁, t₁, t₀)
     retcode = currstatus(tracker)
     if retcode == :success
         x₀ .= currx(tracker)
     end
     retcode
end
function track!(tracker::PathTracker, x₁, t₁, t₀)
    setup!(tracker, x₁, t₁, t₀)

    while currstatus(tracker) == :ok
        step!(tracker)
        check_terminated!(tracker)
    end
    returncode = currstatus(tracker)
    if returncode == :success
        res = refine!(tracker)
    end
    returncode
end

"""
    setup!(pathtracker, x₁, t₁, t₀, checkstartvalue=true))

Setup `pathtracker` to track `x₁` from `t₁` to `t₀`.
"""
function setup!(tracker::PathTracker, x₁, t₁, t₀, checkstartvalue=true)
    S = tracker.state

    S.start = t₁
    S.target = t₀
    StepLength.reset!(S.steplength, tracker.steplength, t₁, t₀)
    S.t = 1.0
    S.Δt = StepLength.relsteplength(S.steplength)
    S.status = :ok

    S.iters = 0
    tracker.x .= x₁

    Homotopies.precondition!(tracker.cache.homotopy, tracker.x, t₁)

    if checkstartvalue
        checkstartvalue!(tracker)
    end
    tracker
end
patch(cache::Cache) = cache.homotopy.homotopy.patch # this is a little bit hacky...

function currresidual(tracker)
    Homotopies.evaluate!(tracker.cache.out, tracker.cache.homotopy, currx(tracker), currt(tracker))
    infinity_norm(tracker.cache.out)
end

function checkstartvalue!(tracker)
    res = currresidual(tracker)
    tol = tracker.options.tol * 20
    if res > tol
        tracker.state.status = :invalid_startvalue
    end
    nothing
end

function step!(tracker, emit_update=true)
    state, cache = tracker.state, tracker.cache

    state.iters += 1
    step_successfull, status = PredictionCorrection.predict_correct!(
        currx(tracker),
        tracker.predictor_corrector,
        cache.predictor_corrector,
        cache.homotopy, currt(state), currΔt(state),
        tracker.options.tol,
        tracker.options.corrector_maxiters)

    update_t_Δt!(tracker.state, tracker.steplength, step_successfull)
    nothing
end

function update_t_Δt!(state, steplength, step_successfull)
    if step_successfull
        state.t -= state.Δt
    end

    status = StepLength.update!(state.steplength, steplength, step_successfull)
    if status != :ok
        state.status = status
        return nothing
    end

    Δt = StepLength.relsteplength(state.steplength)
    # Due to numerical errors we would sometimes be at t=0.1+2eps() and we would just
    # do a step of size 0.1 and then again one of size 2eps() which is quite wasteful.
    # So we correct this here and try to do a minimal larger step.
    if Δt + 5e-16 > state.t
        state.Δt = state.t
    else
        state.Δt = Δt
    end
    nothing
end

function check_terminated!(tracker)
    if tracker.state.t == 0.0
        tracker.state.status = :success
    elseif tracker.state.iters ≥ tracker.options.maxiters
        tracker.state.status = :maxiters
    end
    nothing
end

function refine!(tracker)
    PredictionCorrection.refine!(currx(tracker),
        tracker.predictor_corrector,
        tracker.cache.predictor_corrector,
        tracker.cache.homotopy, currt(tracker),
        tracker.options.refinement_tol,
        tracker.options.refinement_maxiters,
        currresidual(tracker))
end
