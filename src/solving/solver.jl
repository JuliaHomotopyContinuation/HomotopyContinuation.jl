struct SolverOptions
    endgame_start::Float64
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
        E<:Endgaming.Endgame, PS<:Union{Nothing,PatchSwitching.PatchSwitcher}, C<:SolverCache}
    prob::P
    tracker::T
    endgame::E
    patchswitcher::PS
    t₁::Float64
    t₀::Float64
    seed::Int
    options::SolverOptions
    cache::C
end

function Solver(prob::Problems.AbstractProblem, start_solutions, args...; kwargs...) where {T<:AbstractVector}
    Solver(prob, Utilities.start_solution_sample(start_solutions), args...; kwargs...)
end

function Solver(prob::Problems.AbstractProblem, startsolutionsample::AbstractVector{<:Complex}, t₁, t₀, seed=0;
    endgame_start=0.1,
    report_progress=true,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, report_progress)
    Solver(prob, startsolutionsample, t₁, t₀, seed, options; kwargs...)
end

function Solver(prob::Problems.Projective, startsolutionsample::AbstractVector{<:Complex}, t₁, t₀, seed, options::SolverOptions;kwargs...)
    x₁= Problems.embed(prob, startsolutionsample)

    tracker = PathTracking.PathTracker(prob, x₁, t₁, t₀; filterkwargs(kwargs, PathTracking.allowed_keywords)...)

    endgame = Endgaming.Endgame(prob.homotopy, x₁; filterkwargs(kwargs, Endgaming.allowed_keywords)...)

    check_kwargs(kwargs)

    Solver(prob, tracker, endgame, nothing, #patchswitcher(prob, x₁, t₀),
        t₁, t₀, seed, options, SolverCache(prob, tracker))
end

function solver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_keywords)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    Solver(prob, startsolutions, 1.0, 0.0, prob.seed; rest...), startsolutions
end

check_kwargs(kwargs) = check_kwargs_empty(invalid_kwargs(kwargs), allowed_keywords())
allowed_keywords() = [:patch, PathTracking.allowed_keywords..., Endgaming.allowed_keywords...]

function invalid_kwargs(kwargs)
    invalids = []
    allowed = allowed_keywords()
    for kwarg in kwargs
        kw = first(kwarg)
        if !any(isequal(kw), allowed)
            push!(invalids, kwarg)
        end
    end
    invalids
end

function patchswitcher(prob::Problems.Projective{H, Problems.NullHomogenization}, x, t₀) where {H<:Homotopies.AbstractHomotopy}
    nothing
end
function patchswitcher(prob::Problems.Projective, x, t₀)
    p₁ = AffinePatches.state(AffinePatches.FixedPatch(), x)
    p₀ = AffinePatches.state(AffinePatches.EmbeddingPatch(), x)
    PatchSwitching.PatchSwitcher(prob.homotopy, p₁, p₀, x, t₀)
end

const Solvers = Vector{<:Solver}
