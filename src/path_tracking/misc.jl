export currt, currΔt,
    curriters, currstatus,
    currx,
    tol,
    corrector_maxiters,
    refinement_tol,
    refinement_maxiters,
    set_tol!,
    set_corrector_maxiters!,
    set_refinement_tol!,
    set_refinement_maxiters!

function Base.show(io::IO, ::MIME"text/plain", tracker::PathTracker)
    print("PathTracker")
end

"""
    checkstart(H, x)

Check whether the `x` has the correct size.
"""
function checkstart(H, x)
    N = Homotopies.nvariables(H)
    N != length(x) && throw(error("Expected `x` to have length $(N) but `x` has length $(length(x))"))
    nothing
end

"""
     currt(tracker::PathTracker)

Current `t`.
"""
currt(tracker::PathTracker) = currt(tracker.state)
currt(state::State) = (1-state.t) * state.target + state.t * state.start

"""
     Δt(tracker::PathTracker)

Current steplength `Δt`.
"""
currΔt(tracker::PathTracker) = currΔt(tracker.state)
currΔt(state::State) = state.Δt * (state.target - state.start)

"""
     iters(tracker::PathTracker)

Current number of iterations.
"""
curriters(tracker::PathTracker) = curriters(tracker.state)
curriters(state::State) = state.iters

"""
     status(tracker::PathTracker)

Current status.
"""
currstatus(tracker::PathTracker) = currstatus(tracker.state)
currstatus(state::State) = state.status

"""
    currx(tracker::PathTracker)

Return the current value of `x`.
"""
currx(tracker::PathTracker) = tracker.x

"""
     tol(tracker::PathTracker)

Current tolerance.
"""
tol(tracker::PathTracker) = tracker.options.tol

"""
     set_tol!(tracker::PathTracker, tol)

Set the current tolerance to `tol`.
"""
function set_tol!(tracker::PathTracker, tol)
     tracker.options.tol = tol
     tol
end

"""
     refinement_tol(tracker::PathTracker)

Current refinement tolerance.
"""
refinement_tol(tracker::PathTracker) = tracker.options.refinement_tol

"""
     set_refinement_maxiters!(tracker::PathTracker, tol)

Set the current refinement tolerance to `tol`.
"""
function set_refinement_tol!(tracker::PathTracker, tol)
     tracker.options.refinement_tol = tol
     tol
end

"""
     refinement_maxiters(tracker::PathTracker)

Current refinement maxiters.
"""
refinement_maxiters(tracker::PathTracker) = tracker.options.refinement_maxiters

"""
     set_refinement_maxiters!(tracker::PathTracker, n)

Set the current refinement maxiters to `n`.
"""
function set_refinement_maxiters!(tracker::PathTracker, n)
     tracker.options.refinement_maxiters = n
     n
end

"""
     corrector_maxiters(tracker::PathTracker)

Current correction maxiters.
"""
corrector_maxiters(tracker::PathTracker) = tracker.options.corrector_maxiters

"""
     set_corrector_maxiters!(tracker::PathTracker, n)

Set the current correction maxiters to `n`.
"""
function set_corrector_maxiters!(tracker::PathTracker, n)
     tracker.options.corrector_maxiters = n
     n
end
