struct SolverOptions
    endgame_start::Float64
    report_progress::Bool
end

struct SolverCache{P<:PathResultCache}
    pathresult::P
end

function SolverCache(prob, tracker)
    pathresult = PathResultCache(prob, currx(tracker))

    SolverCache(pathresult)
end

struct Solver{P<:AbstractProblem, T<:PathTracker,
        E<:Endgame, C<:SolverCache}
    prob::P
    tracker::T
    endgame::E
    t₁::Float64
    t₀::Float64
    seed::Int
    options::SolverOptions
    cache::C
end

function Solver(prob::AbstractProblem, start_solutions, args...; kwargs...) where {T<:AbstractVector}
    Solver(prob, start_solution_sample(start_solutions), args...; kwargs...)
end

function Solver(prob::AbstractProblem, startsolutionsample::AbstractVector{<:Complex}, t₁, t₀, seed=0;
    endgame_start=0.1,
    report_progress=true,
    kwargs...)
    !(t₀ ≤ endgame_start ≤ t₁) && throw(error("`endgame_start` has to be between `t₁` and`t₀`"))

    options = SolverOptions(endgame_start, report_progress)
    Solver(prob, startsolutionsample, t₁, t₀, seed, options; kwargs...)
end

function Solver(prob::ProjectiveProblem, startsolutionsample::AbstractVector{<:Complex}, t₁, t₀, seed, options::SolverOptions;kwargs...)
    x₁= embed(prob, startsolutionsample)

    tracker = PathTracker(prob, x₁, t₁, t₀; filterkwargs(kwargs, pathtracker_allowed_keywords)...)

    check_at_infinity = homvars(prob) !== nothing
    endgame = Endgame(prob.homotopy, x₁; check_at_infinity=check_at_infinity, filterkwargs(kwargs, endgame_allowed_keywords)...)

    check_kwargs(kwargs)

    Solver(prob, tracker, endgame, t₁, t₀, seed, options, SolverCache(prob, tracker))
end

function solver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    Solver(prob, startsolutions, 1.0, 0.0, prob.seed; rest...), startsolutions
end

check_kwargs(kwargs) = check_kwargs_empty(invalid_kwargs(kwargs), allowed_keywords())
allowed_keywords() = [:patch, pathtracker_allowed_keywords..., endgame_allowed_keywords...]

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
