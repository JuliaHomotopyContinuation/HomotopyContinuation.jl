struct SolverOptions
    endgame_start::Float64
    pathcrossing_tol::Float64
    report_progress::Bool
end

struct SolverCache{P<:PathResultCache}
    pathresult::P
end

function SolverCache(prob, tracker)
    pathresult = PathResultCache(prob, PathTracking.currx(tracker))

    SolverCache(pathresult)
end

struct Solver{P<:Problems.AbstractProblem, T<:PathTracking.PathTracker,
        E<:Endgame.Endgamer, PS<:PatchSwitching.PatchSwitcher, C<:SolverCache}
    prob::P
    tracker::T
    endgamer::E
    patchswitcher::PS
    t₁::Float64
    t₀::Float64
    options::SolverOptions
    cache::C
end


function Solver(prob::Problems.AbstractProblem, start_solutions, t₁, t₀=0.0;
    endgame_start=0.1,
    pathcrossing_tol=1e-6,
    report_progress=true,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, pathcrossing_tol, report_progress)
    Solver(prob, start_solutions, t₁, t₀, options; kwargs...)
end

function Solver(prob::Problems.ProjectiveStartTargetProblem, start_solutions, t₁, t₀, options::SolverOptions; kwargs...)
    x₁ = first(start_solutions)
    @assert x₁ isa AbstractVector

    tracker = pathtracker(prob, x₁, t₁, t₀; kwargs...)
    endgamer = Endgame.Endgamer(tracker, options.endgame_start)
    switcher = patchswitcher(prob, x₁, t₀; kwargs...)

    cache = SolverCache(prob, tracker)
    Solver(prob,
        tracker,
        endgamer,
        switcher,
        t₁, t₀,
        options,
        cache)
end

function pathtracker(prob::Problems.ProjectiveStartTargetProblem, x₁, t₁, t₀; patch=AffinePatches.OrthogonalPatch(), kwargs...)
    x = Problems.embed(prob, x₁)
    H = Homotopies.PatchedHomotopy(prob.homotopy, AffinePatches.state(patch, x))
    PathTracking.PathTracker(H, x, t₁, t₀; kwargs...)
end

function patchswitcher(prob::Problems.ProjectiveStartTargetProblem, x₁, t₀; patch=AffinePatches.OrthogonalPatch(), kwargs...)
    x = Problems.embed(prob, x₁)
    p₁ = AffinePatches.state(patch, x)
    p₀ = AffinePatches.state(AffinePatches.EmbeddingPatch(), x)
    PatchSwitching.PatchSwitcher(prob.homotopy, p₁, p₀, x, t₀)
end

const Solvers = Vector{<:Solver}
