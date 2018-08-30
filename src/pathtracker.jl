export pathtracker_startsolutions

"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a `PathTracking.PathTracker` and `startsolutions` in the same way [`solve`](@ref)
does it. This als takes the same input arguments as `solve`. This is convient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_kwargs)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    t₁, t₀ = one(ComplexF64), zero(ComplexF64)
    embeded_startsolutions = map(x -> Problems.embed(prob, x), startsolutions)
    tracker = PathTracking.PathTracker(prob, embeded_startsolutions[1], t₁, t₀; rest...)

    (tracker=tracker, startsolution=embeded_startsolutions)
end
