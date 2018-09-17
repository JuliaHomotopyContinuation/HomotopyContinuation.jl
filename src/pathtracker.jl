export pathtracker_startsolutions

"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) and `startsolutions` in the same way [`solve`](@ref)
does it. This als0 takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_kwargs)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    t₁, t₀ = one(ComplexF64), zero(ComplexF64)
    tracker = PathTracking.PathTracker(prob, startsolutions[1], t₁, t₀; rest...)

    (tracker=tracker, startsolutions=startsolutions)
end
