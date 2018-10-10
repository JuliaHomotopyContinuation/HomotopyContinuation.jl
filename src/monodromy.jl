module Monodromy

export monodromy_solve

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticArrays: SVector, @SVector
import ..Homotopies, ..PathTracking, ..ProjectiveVectors
using ..Utilities


include("monodromy/group_actions.jl")
include("monodromy/options.jl")
include("monodromy/statistics.jl")
include("monodromy/strategy.jl")

struct MonodromyResult{N, T}
    returncode::Symbol
    solutions::Vector{SVector{N, T}}
    statistics::Statistics
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike},
        p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)

        monodromy_solve(F, p₀, [solution]; kwargs...)
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike},
        p₀::AbstractVector{<:Number},
        solutions::Vector{<:AbstractVector{<:Number}}; kwargs...)
    monodromy_solve(F, SVector{length(p₀)}(p₀), static_solutions(solutions, Val(length(solutions[1]))); kwargs...)
end

function static_solutions(V::Vector{<:AbstractVector{<:Complex{<:AbstractFloat}}}, ::Val{N}) where {N}
    map(v -> SVector{N}(v), V)
end
function static_solutions(V::Vector{<:AbstractVector{<:Number}}, ::Val{N}) where {N}
    map(v -> complex.(float.(SVector{N}(v))), V)
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{T}},
        p₀::SVector{NParams, <:Number},
        solutions::Vector{<:SVector{NVars, <:Complex}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=Triangle(), options...) where {T, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        return error("Number of provided parameters doesn't match the length of initially provided parameter `p₀`.")
    end
    isrealsystem = eltype(p₀) <: Real && T <: Real
    opts = Options(isrealsystem; options...)

    uniquesolutions = UniquePoints(solutions, tol=opts.tol)


    tracker = PathTracking.pathtracker(F, solutions, parameters=parameters, p₁=p₀, p₀=p₀; tol=opts.tol)
    # assemble strategy stuff
    strategy_params, strategy_cache = strategy_parameters_cache(strategy, tracker, p₀)
    statistics = Statistics()

    retcode = monodromy_solve!(
        uniquesolutions,
        statistics,
        tracker,
        p₀,
        strategy_params,
        strategy_cache,
        opts
        )

    MonodromyResult(retcode, points(uniquesolutions), statistics)
end

function strategy_parameters_cache(strategy, tracker, p₀)
    parameters(strategy, p₀), cache(strategy, tracker)
end

function monodromy_solve!(
    solutions::UniquePoints,
    stats::Statistics,
    tracker::PathTracking.PathTracker,
    p₀::SVector,
    parameters::AbstractStrategyParameters,
    strategy_cache::AbstractStrategyCache,
    options::Options)

    t₀ = time_ns()
    k = 0
    n = length(solutions)
    generated_parameters!(stats, n)
    # We prepopulate the solutions
    for i=1:n
        retcode = apply_group_actions_greedily!(solutions, solutions[i], options)
        if retcode == :done
            return :success
        end
    end


    while length(solutions) < options.target_solutions_count
        retcode = track_set!(solutions, stats, tracker, p₀, parameters, strategy_cache, options)

        if retcode == :done
            break
        end

        dt = (time_ns() - t₀) * 1e-9
        if dt > options.timeout
            return :timeout
        end
        parameters = regenerate(parameters)
        generated_parameters!(stats, length(solutions))
    end

    :success
end

function track_set!(solutions::UniquePoints, stats::Statistics, tracker, p₀, params::AbstractStrategyParameters, strategy_cache, options::Options)
    queue = copy(points(solutions))
    while !isempty(queue)
        s₀ = pop!(queue)
        s₁, retcode = loop(tracker, s₀, p₀, params, strategy_cache, stats)
        if retcode == :success && add!(solutions, s₁)
            push!(queue, s₁)
            if options.done_callback(s₁) || length(solutions) ≥ options.target_solutions_count
                return :done
            end

            retcode = apply_group_actions_greedily!(solutions, s₁, options, queue)
            if retcode == :done
                return :done
            end
        end
    end
    :incomplete
end

function apply_group_actions_greedily!(solutions::UniquePoints, s, options, queue=nothing)
    for tᵢ in options.group_actions(s)
        if add!(solutions, tᵢ)
            if queue !== nothing
                push!(queue, tᵢ)
            end
            if options.done_callback(tᵢ) || length(solutions) ≥ options.target_solutions_count
                return :done
            end
            retcode = apply_group_actions_greedily!(solutions, tᵢ, options, queue)
            if retcode == :done
                return :done
            end
        end
    end
    return :incomplete
end

end
