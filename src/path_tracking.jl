module PathTracking

export pathtracker_startsolutions, pathtracker

include("path_tracking/types.jl")
include("path_tracking/misc.jl")
include("path_tracking/tracking.jl")

"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) and `startsolutions` in the same way [`solve`](@ref)
does it. This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_keywords)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    tracker = PathTracking.PathTracker(prob, Utilities.start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); rest...)

    (tracker=tracker, startsolutions=startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) in the same way [`solve`](@ref)
does it. This als0 takes the same input arguments as [`solve`](@ref). This is convenient if you want
to investigate single paths.
"""
function pathtracker(args...; kwargs...)
    tracker, _ = pathtracker_startsolutions(args...; kwargs...)
    tracker
end


end
