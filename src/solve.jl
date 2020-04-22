export solve, solver_startsolutions

struct Solver{T<:AbstractTracker}
    tracker::T
end

function solver_startsolutions(F::System; start_system = :polyhedral, kwargs...)
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

    Solver(tracker), starts
end

"""
    solve(...)

TODO
"""
function solve(args...; kwargs...)
    solver, starts = solver_startsolutions(args...; kwargs...)
    solve(solver, starts)
end
solve(S::Solver, starts) = track.(S.tracker, starts)
