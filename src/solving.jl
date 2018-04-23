module Solving

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ..Endgame
import ..Homotopies
import ..PathTracking
import ..Problems
import ..Systems

export Solver,
    solve


struct SolverOptions
    endgame_start::Float64
    pathcrossing_tol::Float64
end


struct Solver{P<:Problems.AbstractProblem, T<:PathTracking.PathTracker, E<:Endgame.Endgamer, V<:AbstractVector}
    prob::P
    tracker::T
    endgamer::E
    t₁::Float64
    t₀::Float64
    start_solutions::Vector{V}
    options::SolverOptions
end


function Solver(prob::Problems.AbstractProblem, start_solutions, t₁, t₀=0.0;
    endgame_start=0.1,
    pathcrossing_tol=1e-6,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, pathcrossing_tol)
    Solver(prob, start_solutions, t₁, t₀, options; kwargs...)
end

function Solver(prob::Problems.ProjectiveStartTargetProblem, start_solutions, t₁, t₀, options::SolverOptions; kwargs...)
    @assert !isempty(start_solutions) "`start_solutions` are empty"
    if start_solutions isa Vector
        x₁s = start_solutions
    else
        x₁s = collect(start_solutions)
    end
    x₀ = first(x₁s)
    @assert x₀ isa AbstractVector

    tracker = PathTracking.PathTracker(prob, x₀, t₁, t₀; kwargs...)
    endgamer = Endgame.Endgamer(tracker, options.endgame_start)
    Solver(prob,
        tracker,
        endgamer,
        t₁, t₀,
        x₁s,
        options)
end

include("solving/path_crossing.jl")
include("solving/path_result.jl")

function solve(solver::Solver)
    # We set the number of BLAS threads to 1 since we multithread by ourself.
    # But even if we would not, the overhead of the threading is so large
    # that the single threaded version is around two times faster (at least on Mac Julia v0.6.2)
    nblas_threads = Utilities.get_num_BLAS_threads()
    BLAS.set_num_threads(1)

    t_endgame = solver.options.endgame_start
    # Phase 1 Track until endgame zone
    pre_endgame_paths = dmap(solver.options, solver.start_solutions) do x₀
        trackpath(solver, x₀, solver.t₁, t_endgame)
    end

    if solver.t₀ == t_endgame
        BLAS.set_num_threads(nblas_threads)
        # return pre_endgame_paths
        return pathresults(solver, pre_endgame_paths)
    end

    ncrossedpaths, crossing_indices = pathcrossing_check!(pre_endgame_paths, solver)

    BLAS.set_num_threads(nblas_threads)

    endgame_results = dmap(solver.options, pre_endgame_paths) do r
        Endgame.play(solver.endgamer, r.x.data, t_endgame)
    end
    # pre_endgame_paths
    return pathresults(solver, endgame_results)
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
        dest[i] = f(src[i])
    end
    dest
end




end
