module Monodromy

export monodromy_solve

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticArrays: SVector, @SVector
import ..Homotopies, ..PathTracking, ..ProjectiveVectors


mutable struct MonodromyStatistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    ninstances::Int
end
MonodromyStatistics() = MonodromyStatistics(0, 0, 0)

function pathtracked!(stats::MonodromyStatistics, retcode)
    if retcode == :success
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end

struct MonodromyResult{T}
    returncode::Symbol
    solutions::Vector{Vector{T}}
    statistics::MonodromyStatistics
end


############
# Strategy
############

include("monodromy/strategy.jl")
include("monodromy/group_actions.jl")


##############
# Statistics
##############




struct MonodromyOptions{F1<:Function, F2<:Tuple}
    target_solutions_count::Int
    timeout::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
end

function MonodromyOptions(;
    target_solutions_count=error("target_solutions_count not provided"),
    timeout=float(typemax(Int)),
    done_callback=always_false,
    group_action=nothing,
    group_actions=GroupActions(group_action))

    MonodromyOptions(target_solutions_count, float(timeout), done_callback, GroupActions(group_actions))
end

always_false(x) = false
complex_conjugation(x) = (conj.(x),)




function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike},
        p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)

        monodromy_solve(F, p₀, [solution]; kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{T}},
        p₀::AbstractVector{<:Number},
        solutions::Vector{<:Vector{<:Number}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=Triangle(), options...) where {T}

    comp_solutions = map(v -> complex.(v), solutions)
    tracker = PathTracking.pathtracker(F, comp_solutions, parameters=parameters, p₁=p₀, p₀=p₀)
    # assemble strategy stuff
    nparams = length(p₀)
    p₀ = SVector{nparams}(p₀)
    strategy_params, strategy_cache = strategy_parameters_cache(strategy, tracker, p₀)
    statistics = MonodromyStatistics()
    retcode = monodromy_solve!(
        comp_solutions,
        statistics,
        tracker,
        p₀,
        strategy_params,
        strategy_cache,
        MonodromyOptions(;options...)
        )

    MonodromyResult(retcode, comp_solutions, statistics)
end

function strategy_parameters_cache(strategy, tracker, p₀)
    parameters(strategy, p₀), cache(strategy, tracker)
end

function monodromy_solve!(
    solutions::Vector{<:Vector{<:Complex}},
    stats::MonodromyStatistics,
    tracker::PathTracking.PathTracker,
    p₀::SVector,
    parameters::AbstractStrategyParameters,
    strategy_cache::AbstractStrategyCache,
    options::MonodromyOptions)

    t₀ = time_ns()
    k = 0

    # We prepopulate the solutions
    n = length(solutions)
    for i=1:n
        retcode = apply_group_actions_greedily!(nothing, solutions, solutions[i], options)
        if retcode == :done
            return :success
        end
    end


    while length(solutions) < options.target_solutions_count
        stats.ninstances += 1
        retcode = track_set!(solutions, stats, tracker, p₀, parameters, strategy_cache, options)

        if retcode == :done
            break
        end

        dt = (time_ns() - t₀) * 1e-9
        if dt > options.timeout
            return :timeout
        end
        parameters = regenerate(parameters)
    end

    :success
end

function track_set!(solutions, stats::MonodromyStatistics, tracker, p₀, params::AbstractStrategyParameters, strategy_cache, options::MonodromyOptions)
    S = copy(solutions)
    while !isempty(S)
        s₀ = pop!(S)
        s₁, retcode = loop(tracker, s₀, p₀, params, strategy_cache, stats)
        if retcode == :success && !iscontained(solutions, s₁)
            push!(solutions, s₁)
            push!(S, s₁)
            if options.done_callback(s₁) || length(solutions) ≥ options.target_solutions_count
                return :done
            end

            retcode = apply_group_actions_greedily!(S, solutions, s₁, options)
            if retcode == :done
                return :done
            end
        end
    end
    :incomplete
end

function apply_group_actions_greedily!(S, solutions, s, options)
    for tᵢ in options.group_actions(s)
        if !iscontained(solutions, tᵢ)
            push!(solutions, tᵢ)
            if !(S isa Nothing)
                push!(S, tᵢ)
            end
            if options.done_callback(tᵢ) || length(solutions) ≥ options.target_solutions_count
                return :done
            end
            retcode = apply_group_actions_greedily!(S, solutions, tᵢ, options)
            if retcode == :done
                return :done
            end
        end
    end
    return :incomplete
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

end
