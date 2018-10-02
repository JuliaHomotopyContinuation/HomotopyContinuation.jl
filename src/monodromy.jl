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


struct MonodromyOptions{F<:Function}
    target_solutions_count::Int
    timeout::Float64
    done_callback::F
end

function MonodromyOptions(;
    target_solutions_count=error("target_solutions_count not provided"),
    timeout=float(typemax(Int)),
    done_callback=always_false)

    MonodromyOptions(target_solutions_count, float(timeout), done_callback)
end


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
        strategy=Triangle(), options...) where {T}

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
        MonodromyOptions(;options...)
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
    options::MonodromyOptions)

    t₀ = time_ns()
    k = 0
    while length(solutions) < options.target_solutions_count
        retcode = track_set!(solutions, tracker, p₀, parameters, strategy_cache, options)

        if retcode == :done
            break
        end

        dt = (time_ns() - t₀) * 1e-9
        if dt > options.timeout
            break
        end
        parameters = regenerate(parameters)
    end
    solutions
end

function track_set!(solutions, tracker, p₀, params::MonodromyStrategyParameters, strategy_cache, options::MonodromyOptions)
    S = copy(solutions)
    while !isempty(S)
        s₀ = pop!(S)
        s₁, retcode = loop(tracker, s₀, p₀, params, strategy_cache)
        if retcode == :success && !iscontained(solutions, s₁)
            push!(solutions, s₁)
            push!(S, s₁)
            if options.done_callback(s₁) || length(solutions) ≥ options.target_solutions_count
                return :done
            end
        end
    end
    :incomplete
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

always_false(x) = false
