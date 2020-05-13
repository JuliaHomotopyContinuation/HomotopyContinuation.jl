export solve, Solver, solver_startsolutions

struct Solver{T<:AbstractPathTracker}
    tracker::T
    seed::UInt32
end
Solver(tracker) = Solver(tracker, UInt32(0)) 

function solver_startsolutions(
    F::Union{System,AbstractSystem};
    seed = rand(UInt32),
    start_system = :polyhedral,
    kwargs...,
)
    Random.seed!(seed)
    if start_system == :polyhedral
        tracker, starts = polyhedral(F; kwargs...)
    elseif start_system == :total_degree
        tracker, starts = total_degree(F; kwargs...)
    else
        throw(KeywordArgumentException(
            :start_system,
            start_system,
            "Possible values are: `:polyhedral` and `:total_degree`.",
        ))
    end

    Solver(tracker, seed), starts
end

"""
    solve(...)

TODO
"""
function solve(args...; kwargs...)
    solver, starts = solver_startsolutions(args...; kwargs...)
    solve(solver, starts)
end
function solve(S::Solver, starts)
    path_results = PathResult[]
    for (k, s) in enumerate(starts)
        push!(path_results, track(S.tracker, s; path_number = k))
    end

    Result(path_results; seed = S.seed)
end
