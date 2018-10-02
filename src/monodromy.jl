export monodromy_solve

import StaticArrays: SVector, @SVector
import LinearAlgebra


############
# Strategy
############

include("monodromy/strategy.jl")


##############
# Statistics
##############

struct MonodromyStatistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    ntriangles::Int
end
MonodromyStatistics() = MonodromyStatistics(0, 0, 0)

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike},
        p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)

        monodromy_solve(F, p₀, [solution]; kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{T}},
        p₀::AbstractVector{<:Number},
        solutions::Vector{<:Vector{<:Number}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        target_solutions_count=error("target_solutions_count not provided"),
        strategy=Triangle()) where {T}

    comp_solutions = map(v -> complex.(v), solutions)
    tracker = pathtracker(F, comp_solutions, parameters=parameters, p₁=p₀, p₀=p₀)
    # assemble strategy stuff
    nparams = length(p₀)
    p₀ = SVector{nparams}(p₀)
    strategy_params, strategy_cache = strategy_parameters_cache(strategy, tracker, p₀)
    monodromy_solve!(
        comp_solutions,
        tracker,
        p₀,
        strategy_params,
        strategy_cache,
        target_solutions_count
        )
end

function strategy_parameters_cache(strategy, tracker, p₀)
    parameters(strategy, p₀), cache(strategy, tracker)
end

function monodromy_solve!(
    solutions::Vector{<:Vector{<:Complex}},
    tracker::PathTracking.PathTracker,
    p₀::SVector,
    parameters::MonodromyStrategyParameters,
    strategy_cache::MonodromyStrategyCache,
    target_solutions_count::Integer)

    t₀ = time_ns()
    k = 0
    while length(solutions) < target_solutions_count
        track_set!(solutions, tracker, p₀, parameters, strategy_cache)

        dt = (time_ns() - t₀) * 1e-9
        if dt > 10
            break
        end
        parameters = regenerate(parameters)
    end
    solutions
end

function track_set!(solutions, tracker, p₀, params::MonodromyStrategyParameters, strategy_cache)
    S = copy(solutions)
    while !isempty(S)
        s₀ = pop!(S)
        s₁, retcode = loop(tracker, s₀, p₀, params, strategy_cache)
        if retcode == :success && !iscontained(solutions, s₁)
            push!(solutions, s₁)
            push!(S, s₁)
        end
    end
    solutions
end

function iscontained(solutions::Vector{T}, s_new; tol=1e-5) where {T<:AbstractVector}
    for s in solutions
        if LinearAlgebra.norm(s - s_new) < tol
        # if Distances.evaluate(Euclidean(), s, s_new) < tol
            return true
        end
    end
    false
end
