module Monodromy

export monodromy_solve

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticArrays: SVector, @SVector
import ..Homotopies, ..PathTracking, ..ProjectiveVectors, ..Utilities

import ..Utilities: UniquePoints


include("monodromy/group_actions.jl")
include("monodromy/options.jl")
include("monodromy/statistics.jl")
include("monodromy/graph.jl")
include("monodromy/strategy.jl")
include("monodromy/jobs.jl")

struct MonodromyResult{N, T}
    returncode::Symbol
    solutions::Vector{SVector{N, T}}
    statistics::Statistics
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)
    monodromy_solve(F, p₀, [solution]; kwargs...)
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, p₀::AbstractVector{<:Number}, solutions::Vector{<:AbstractVector{<:Number}}; kwargs...)
    monodromy_solve(F, SVector{length(p₀)}(p₀), static_solutions(solutions); kwargs...)
end


static_solutions(V::Vector) = static_solutions(V, Val(length(V[1])))
function static_solutions(V::Vector, ::Val{N}) where {N}
    map(v -> complex.(float.(SVector{N}(v))), V)
end
function static_solutions(V::Vector{<:AbstractVector{<:Complex{<:AbstractFloat}}}, ::Val{N}) where {N}
    SVector{N}.(V)
end

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{T}},
        p₀::SVector{NParams, <:Number},
        startsolutions::Vector{<:SVector{NVars, <:Complex}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=Triangle(), options...) where {T, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        return error("Number of provided parameters doesn't match the length of initially provided parameter `p₀`.")
    end
    isrealsystem = eltype(p₀) <: Real && T <: Real
    opts = Options(isrealsystem; options...)

    #assemble
    G = graph(strategy, p₀, first(startsolutions), opts)
    add_initial_solutions!(G, startsolutions; tol=opts.tol)
    tracker = PathTracking.pathtracker(F, startsolutions, parameters=parameters, p₁=p₀, p₀=p₀; tol=opts.tol)
    statistics = Statistics()
    generated_parameters!(statistics, length(solutions(G)))
    # solve
    retcode = monodromy_solve!(G, tracker, opts, statistics)

    MonodromyResult(retcode, Utilities.points(solutions(G)), statistics)
end

function strategy_parameters_cache(strategy, tracker, p₀)
    parameters(strategy, p₀), cache(strategy, tracker)
end

function monodromy_solve!(G::Graph,
    tracker::PathTracking.PathTracker, options::Options,
    stats::Statistics)

    t₀ = time_ns()

    # # We prepopulate the solutions
    # for i=1:n
    #     retcode = apply_group_actions_greedily!(solutions, solutions[i], options)
    #     if retcode == :done
    #         return :success
    #     end
    # end

    queue = map(x -> Job(x, G.loop[1]), solutions(G))
    iterations_without_progress = 0 # stopping heuristic
    n = length(solutions(G))
    while n < options.target_solutions_count
        retcode = empty_queue!(queue, G, tracker, options, t₀, stats)

        if retcode == :done
            break
        elseif retcode == :timeout
            return :timeout
        end

        n_new = length(solutions(G))
        if n == n_new
            iterations_without_progress += 1
        else
            n = n_new
        end
        if iterations_without_progress == options.maximal_number_of_iterations_without_progress
            return :heuristic_stop
        end

        regenerate!(queue, G, options, stats)
    end

    :success
end



function empty_queue!(queue, G::Graph, tracker::PathTracking.PathTracker, options::Options, t₀::UInt, stats::Statistics)
    while !isempty(queue)
        job = pop!(queue)
        job_result = execute(job, G, tracker, stats)
        if process!(queue, job, job_result, G, options) == :done
            return :done
        end

        # check timeout
        if ns_to_s(time_ns() - t₀) > options.timeout
            return :timeout
        end
    end
    :incomplete
end

ns_to_s(s) = s * 1e-9

function regenerate!(queue, G::Graph, options::Options, stats::Statistics)
    sols = solutions(G)

    # create a new graph by regenerating the parameters (but don't touch our
    # main node)
    regenerate!(G, options.parameter_sampler, stats)


    for x in sols
        push!(queue, Job(x, G.loop[1]))
    end
    nothing
end

end
