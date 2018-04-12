module Solving

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ..NewHomotopies
import ..PathTracking
import ..Problems
import ..Systems

export Solver,
    solve


struct SolverOptions
    endgame_start::Float64
    pathcrossing_tol::Float64
end


struct Solver{P<:Problems.AbstractProblem, M, S, C, SS}
    prob::P
    tracker::PathTracking.PathTracker{M, S, C}
    #endgamer::
    t₁::Float64
    t₀::Float64
    start_solutions::SS
    options::SolverOptions
end


function Solver(prob::Problems.AbstractProblem, start_solutions, t₁, t₀;
    endgame_start=0.1,
    pathcrossing_tol=1e-10,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, pathcrossing_tol)
    Solver(prob, start_solutions, t₁, t₀, options; kwargs...)
end

function Solver(prob::Problems.ProjectiveStartTargetProblem, start_solutions, t₁, t₀, options::SolverOptions; kwargs...)
    @assert !isempty(start_solutions) "`start_solutions` are empty"
    x₀ = first(start_solutions)
    @assert x₀ isa AbstractVector

    tracker = PathTracking.PathTracker(prob, x₀, t₁, t₀; kwargs...)
    Solver(prob,
        tracker,
        t₁, t₀,
        start_solutions,
        options)
end

include("solving/path_crossing.jl")

function solve(solver::Solver)
    t_endgame = solver.options.endgame_start
    # Phase 1 Track until endgame zone
    pre_endgame_paths = dmap(solver.options, solver.start_solutions) do x₀
        trackpath(solver, x₀, solver.t₁, t_endgame)
    end

    if solver.t₀ == t_endgame
        return pre_endgame_paths
    end

    ncrossedpaths, crossing_indices = pathcrossing_check!(pre_endgame_paths, solver)

    pre_endgame_paths
end


trackpath(solver, x₀, t₁, t₀) = PathTracking.track(solver.tracker, Problems.embed(solver.prob, x₀), t₁, t₀)

"""
    dmap(f, options::SolverOptions, c)

Map over collection `c` and apply `f` on each element. This uses the parallelization
strategy defined in  `options`.
"""
function dmap(f::F, options::SolverOptions, src) where {F<:Function}
    #TODO
    map(f, src)
end

"""
    dmap_indices!(f, dest, options::SolverOptions, indices, c)

Map over indices `indices`, apply `f` on each `c[i]` where `i ∈ indices`
and store result in `dest[i]`.
This uses the parallelization strategy defined in  `options`.
"""
function dmap_indices!(f::F, dest, options::SolverOptions, indices, src) where {F<:Function}
    for i in indices
        dest[i] = f([i])
    end
    dest
end




end
