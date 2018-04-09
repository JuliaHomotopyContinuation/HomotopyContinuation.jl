module PathTrackers

import ..NewHomotopies
import ..NewHomotopies: AbstractHomotopy, AbstractStartTargetHomotopy, AbstractHomotopyCache, HomotopyWithCache
import ..Predictors
import ..Correctors
import ..AffinePatches
import ..AffinePatches: AbstractAffinePatch
import ..StepSize
import ..StepSize: AbstractStepSize, AbstractStepSizeState
import ..Problems: ProjectiveStartTargetProblem
import ..Problems

export AbstractPathTracker,
    AbstractPathTrackerCache,
    AbstractPathTrackerState,
    Options

abstract type AbstractPathTracker{H<:AbstractHomotopy} end
abstract type AbstractPathTrackerCache end
abstract type AbstractPathTrackerState end


include("path_trackers/interface.jl")
include("path_trackers/patched_homotopy.jl")
include("path_trackers/predictor_corrector.jl")
include("path_trackers/projective_tracker.jl")


function pathtracker(H::AbstractHomotopy, x, start, target; tracker=ProjectiveTracker, kwargs...)
    pathtracker(tracker(H; kwargs...), x, start, target)
end

function pathtracker(tracker::AbstractPathTracker{<:AbstractHomotopy{M, N}}, x, start, target) where {M, N}
    @assert length(x) == N "Expected initial solution to have length $N but got $(length(x))."
    tracker_state = state(tracker, x, start, target)
    tracker_cache = cache(tracker, tracker_state)
    TrackerIterator(tracker, tracker_state, tracker_cache)
end


end
