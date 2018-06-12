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
        E<:Endgaming.Endgame, PS<:PatchSwitching.PatchSwitcher, C<:SolverCache}
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

function Solver(prob::Problems.AbstractProblem, start_solutions, t₁, t₀, seed=0;
    endgame_start=0.1,
    report_progress=true,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, report_progress)
    Solver(prob, start_solutions, t₁, t₀, seed, options; kwargs...)
end

function Solver(prob::Problems.ProjectiveStartTargetProblem, start_solutions, t₁, t₀, seed, options::SolverOptions; kwargs...)
    x₁ = first(start_solutions)
    @assert x₁ isa AbstractVector
    x = Problems.embed(prob, x₁)

    check_kwargs(kwargs)

    tracker = pathtracker(prob, x, t₁, t₀; kwargs...)
    endgame = Endgaming.Endgame(prob.homotopy, x; kwargs...)
    switcher = patchswitcher(prob, x, t₀)

    cache = SolverCache(prob, tracker)
    Solver(prob,
        tracker,
        endgame,
        switcher,
        t₁, t₀,
        seed,
        options,
        cache)
end

allowed_kwargs() = [:patch, PathTracking.allowed_kwargs...,
    Endgaming.allowed_kwargs...]

function check_kwargs(kwargs)
    invalids = invalid_kwargs(kwargs)
    if !isempty(invalids)
        msg = "Unexpected keyword argument(s): "
        first_el = true
        for kwarg in invalids
            if !first_el
                msg *= ", "
            end
            msg *= "$(first(kwarg))=$(last(kwarg))"
            first_el = false
        end
        msg *= "\nAllowed keywords are\n"
        msg *= join(allowed_kwargs(), ", ")
        throw(ErrorException(msg))
    end
end
function invalid_kwargs(kwargs)
    invalids = []
    allowed = allowed_kwargs()
    for kwarg in kwargs
        kw = first(kwarg)
        if !any(equalto(kw), allowed)
            push!(invalids, kwarg)
        end
    end
    invalids
end

function pathtracker(prob::Problems.ProjectiveStartTargetProblem, x, t₁, t₀; patch=AffinePatches.OrthogonalPatch(), kwargs...)
    H = Homotopies.PatchedHomotopy(prob.homotopy, AffinePatches.state(patch, x))
    PathTracking.PathTracker(H, x, complex(t₁), complex(t₀); kwargs...)
        # filterkwargs(kwargs, PathTracking.PATH_TRACKER_KWARGS)...)
end

function patchswitcher(prob::Problems.ProjectiveStartTargetProblem, x, t₀)
    p₁ = AffinePatches.state(AffinePatches.FixedPatch(), x)
    p₀ = AffinePatches.state(AffinePatches.EmbeddingPatch(), x)
    PatchSwitching.PatchSwitcher(prob.homotopy, p₁, p₀, x, t₀)
end

const Solvers = Vector{<:Solver}
