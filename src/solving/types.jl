struct SolverOptions
    endgame_start::Float64
    pathcrossing_tol::Float64
end

struct SolverCache{P<:PathResultCache}
    pathresult::P
end

function SolverCache(prob, tracker)
    pathresult = PathResultCache(prob, PathTracking.currx(tracker))

    SolverCache(pathresult)
end


struct Solver{P<:Problems.AbstractProblem, T<:PathTracking.PathTracker, E<:Endgame.Endgamer, C<:SolverCache}
    prob::P
    tracker::T
    endgamer::E
    t₁::Float64
    t₀::Float64
    options::SolverOptions
    cache::C
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
    x₀ = first(start_solutions)
    @assert x₀ isa AbstractVector

    tracker = PathTracking.PathTracker(prob, x₀, t₁, t₀; kwargs...)
    endgamer = Endgame.Endgamer(tracker, options.endgame_start)
    cache = SolverCache(prob, tracker)
    Solver(prob,
        tracker,
        endgamer,
        t₁, t₀,
        options,
        cache)
end

const Solvers = Vector{<:Solver}
