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
include("monodromy/loop.jl")

struct MonodromyResult{N, T}
    returncode::Symbol
    solutions::Vector{SVector{N, T}}
    statistics::Statistics
end

########
# Setup
########
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)
    monodromy_solve(F, p₀, [solution]; kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, p₀::AbstractVector{<:Number}, solutions::Vector{<:AbstractVector{<:Number}}; kwargs...)
    monodromy_solve(F, SVector{length(p₀)}(p₀), static_solutions(solutions); kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{TC}},
        p₀::SVector{NParams, TP},
        startsolutions::Vector{<:SVector{NVars, <:Complex}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=default_strategy(TC, TP),
        optionskwargs...) where {TC, TP, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        error("Number of provided parameters doesn't match the length of initially provided parameter `p₀`.")
    end
    options = begin
        isrealsystem = TC <: Real && TP <: Real
        Options(isrealsystem; optionskwargs...)
    end

    #assemble
    loop = Loop(strategy, p₀, startsolutions, options)
    tracker = PathTracking.pathtracker(
        F, startsolutions; parameters=parameters, p₁=p₀, p₀=p₀, tol=options.tol
    )
    statistics = Statistics(nsolutions(loop))

    # solve
    retcode = :not_assigned
    try
        retcode = monodromy_solve!(loop, tracker, options, statistics)
    catch e
        if (e isa InterruptException)
            println("interrupt")
            retcode = :interrupt
        else
            rethrow(e)
        end
    end

    MonodromyResult(retcode, Utilities.points(solutions(loop)), statistics)
end

default_strategy(coeff::Type{<:Number}, p::Type{<:Real}) = Triangle(useweights=true)
default_strategy(coeff::Type{<:Number}, p::Type{<:Number}) = Triangle(useweights=false)

# convert vector of vectors to vector of svectors
static_solutions(V::Vector) = static_solutions(V, Val(length(V[1])))
function static_solutions(V::Vector, ::Val{N}) where {N}
    map(v -> complex.(float.(SVector{N}(v))), V)
end
function static_solutions(V::Vector{<:AbstractVector{<:Complex{<:AbstractFloat}}}, ::Val{N}) where {N}
    SVector{N}.(V)
end


##############
# Actual work
##############

"""
    Job{N, T}

A `Job` is consisting of an `Edge` and a solution to the start node of this edge.
"""
struct Job{N, T}
    x::SVector{N, T}
    edge::Edge
end


function monodromy_solve!(loop::Loop,
    tracker::PathTracking.PathTracker, options::Options,
    stats::Statistics)

    t₀ = time_ns()
    iterations_without_progress = 0 # stopping heuristic
    # intialize job queue
    queue = map(x -> Job(x, loop.edges[1]), solutions(loop))

    n = nsolutions(loop)
    while n < options.target_solutions_count
        retcode = empty_queue!(queue, loop, tracker, options, t₀, stats)

        if retcode == :done
            break
        elseif retcode == :timeout
            return :timeout
        end

        # Iterations heuristic
        n_new = nsolutions(loop)
        n == n_new && (iterations_without_progress += 1)
        n = n_new
        if iterations_without_progress == options.maximal_number_of_iterations_without_progress
            return :heuristic_stop
        end

        regenerate_loop_and_schedule_jobs!(queue, loop, options, stats)
    end

    :success
end

function empty_queue!(queue, loop::Loop, tracker::PathTracking.PathTracker, options::Options, t₀::UInt, stats::Statistics)
    while !isempty(queue)
        job = pop!(queue)
        if process!(queue, job, tracker, loop, options, stats) == :done
            return :done
        end

        # check timeout
        if (time_ns() - t₀) > options.timeout * 1e9 # convert s to ns
            return :timeout
        end
    end
    :incomplete
end

function process!(queue::Vector{<:Job}, job::Job, tracker, loop::Loop, options::Options, stats::Statistics)
    y, retcode = track(tracker, job.x, job.edge, loop, stats)
    if retcode ≠ :success
        return :incomplete
    end
    node = loop.nodes[job.edge.target]
    if !iscontained(node, y, tol=options.tol)
        unsafe_add!(node, y)
        # Check if we are done
        if isdone(node, y, options)
            return :done
        end
        next_edge = nextedge(loop, job.edge)
        push!(queue, Job(y, next_edge))

        # Handle group actions
        # Things are setup up such that for nodes where we want to apply
        # group acitons `node.points !== nothing`
        if node.points !== nothing
            for yᵢ in options.group_actions(y)
                if !iscontained(node, yᵢ, tol=options.tol)
                    unsafe_add!(node, yᵢ)
                    # Check if we are done
                    if isdone(node, yᵢ, options)
                        return :done
                    end
                    push!(queue, Job(yᵢ, next_edge))
                end
            end
        end
    end
    return :incomplete
end

function isdone(node::Node, x, options::Options)
    !node.main_node && return false

    options.done_callback(x) ||
    length(node.points) ≥ options.target_solutions_count
end

function regenerate_loop_and_schedule_jobs!(queue, loop::Loop, options::Options, stats::Statistics)
    sols = solutions(loop)
    # create a new loop by regenerating the parameters (but don't touch our
    # main node)
    regenerate!(loop, options, stats)
    for x in sols
        push!(queue, Job(x, loop.edges[1]))
    end
    nothing
end


end
