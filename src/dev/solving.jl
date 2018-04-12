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


function Solver(prob::Problems.AbstractProblem, start_solutions, t₁, t₀; endgame_start=0.1, kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error(`endgame_start` has to be between `t₁` and`t₀`))

    options = SolverOptions(endgame_start)
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

function solve(solver::Solver)
    t_endgame = solver.options.endgame_start
    # Phase 1 Track until endgame zone
    map_distributed(solver.options, solver.start_solutions) do x₀
        PathTracking.track(solver.tracker, Problems.embed(solver.prob, x₀), solver.t₁, t_endgame)
    end

    
end


"""
    map_distributed(f, options::SolverOptions, c)

Map over collection `c` and apply `f` on each element. This uses the parallelization
strategy defined in  `solveroptions`.
"""
function map_distributed(f::F, options::SolverOptions, src) where {F<:Function}
    #TODO
    map(f, src)
end





end
