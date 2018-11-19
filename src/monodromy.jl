module Monodromy

export monodromy_solve

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import ProgressMeter
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

Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MonodromyResult) = x
function Base.show(io::IO, result::MonodromyResult{N, T}) where {N, T}
    println(io, "MonodromyResult{$N, $T} with $(length(result.solutions)) solutions:")
    println(io, "returncode → $(result.returncode)")
    print(io, "statistics → $(result.statistics)")
end


########
# Setup
########
"""
    monodromy_solve(F, sols, p; parameters=..., options..., pathtrackerkwargs...)

Solve a polynomial system `F(x;p)` with specified parameters and initial solutions `sols`
by monodromy techniques. This makes loops in the parameter space of `F` to find new solutions.

## Options
* `target_solutions_count=nothing`: The computations are stopped if this number of solutions is reached.
* `done_callback=always_false`: A callback to end the computation early. This function takes 2 arguments. The first one is the new solution `x` and the second one are all current solutions (including `x`). Return `true` if the compuation is done.
* `maximal_number_of_iterations_without_progress::Int=10`: The maximal number of iterations (i.e. loops generated) without any progress.
* `group_action=nothing`: A function taking one solution and returning other solutions if there is a constructive way to obtain them, e.g. by symmetry.
* `strategy`: The strategy used to create loops. By default this will be `Trianlge` with weights if `F` is a real system.
* `showprogress=true`: Enable a progress meter.
* `tol::Float64=1e-6`: The tolerance with which paths are tracked and with which it is decided whether to solutions are identical.
* `group_actions=GroupActions(group_action)`: If there is more than one group action you can use this to chain the application of them.
* `group_action_on_all_nodes=false`: By default the group_action(s) are only applied on the solutions with the main parameter `p`. If this is enabled then it is applied for every parameter `q`.
* `parameter_sampler=independent_normal`: A function taking the parameter `p` and returning a new random parameter `q`. By default each entry of the parameter vector is drawn independently from the unviraite normal distribution.
* `timeout=float(typemax(Int))`: The maximal number of *seconds* the computation is allowed to run.
* `minimal_number_of_solutions`: The minimal number of solutions before a stopping heuristic is applied. By default this is half of `target_solutions_count` if applicable otherwise 2.
"""
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, solution::Vector{<:Number}, p₀::AbstractVector{<:Number}; kwargs...)
    monodromy_solve(F, [solution], p₀; kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike}, solutions::Vector{<:AbstractVector{<:Number}}, p₀::AbstractVector{<:Number}; kwargs...)
    monodromy_solve(F, static_solutions(solutions), SVector{length(p₀)}(p₀); kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{TC}},
        startsolutions::Vector{<:SVector{NVars, <:Complex}},
        p₀::SVector{NParams, TP};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=default_strategy(TC, TP),
        showprogress=true,
        kwargs...) where {TC, TP, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        error("Number of provided parameters doesn't match the length of initially provided parameter `p₀`.")
    end

    optionskwargs, restkwargs = Utilities.splitkwargs(kwargs, options_allowed_keywords)
    options = begin
        isrealsystem = TC <: Real && TP <: Real
        Options(isrealsystem; optionskwargs...)
    end

    #assemble
    loop = Loop(strategy, p₀, startsolutions, options)
    tracker = PathTracking.pathtracker(
        F, startsolutions; parameters=parameters, p₁=p₀, p₀=p₀,
        restkwargs...)
    statistics = Statistics(nsolutions(loop))

    # solve
    retcode = :not_assigned
    if showprogress
        progress = ProgressMeter.ProgressUnknown("Solutions found:")
    else
        progress = nothing
    end
    try
        retcode = monodromy_solve!(loop, tracker, options, statistics, progress)
    catch e
        if (e isa InterruptException)
            println("interrupt")
            retcode = :interrupt
        else
            rethrow(e)
        end
    end
    finished!(statistics, nsolutions(loop))
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
    stats::Statistics, progress)

    t₀ = time_ns()
    iterations_without_progress = 0 # stopping heuristic
    # intialize job queue
    queue = map(x -> Job(x, loop.edges[1]), solutions(loop))

    n = nsolutions(loop)
    while n < options.target_solutions_count
        retcode = empty_queue!(queue, loop, tracker, options, t₀, stats, progress)

        if retcode == :done
            update_progress!(progress, loop, stats; finish=true)
            break
        elseif retcode == :timeout
            return :timeout
        end

        # Iterations heuristic
        n_new = nsolutions(loop)
        if n == n_new
            iterations_without_progress += 1
        else
            iterations_without_progress = 0
            n = n_new
        end
        if iterations_without_progress == options.maximal_number_of_iterations_without_progress &&
            n_new ≥ options.minimal_number_of_solutions
            return :heuristic_stop
        end

        regenerate_loop_and_schedule_jobs!(queue, loop, options, stats)
    end

    :success
end

function empty_queue!(queue, loop::Loop, tracker::PathTracking.PathTracker, options::Options,
        t₀::UInt, stats::Statistics, progress)
    while !isempty(queue)
        job = pop!(queue)
        if process!(queue, job, tracker, loop, options, stats, progress) == :done
            return :done
        end
        update_progress!(progress, loop, stats)
        # check timeout
        if (time_ns() - t₀) > options.timeout * 1e9 # convert s to ns
            return :timeout
        end
    end
    :incomplete
end

function process!(queue::Vector{<:Job}, job::Job, tracker, loop::Loop, options::Options,
                  stats::Statistics, progress)
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
        # group actions `node.points !== nothing`
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

function update_progress!(::Nothing, loop::Loop, statistics::Statistics; finish=false)
    nothing
end
function update_progress!(progress, loop::Loop, statistics::Statistics; finish=false)
    ProgressMeter.update!(progress, length(solutions(loop)), showvalues=(
        ("# paths tracked: ", statistics.ntrackedpaths),
        ("# loops generated: ", statistics.nparametergenerations)
    ))
    if finish
        ProgressMeter.finish!(progress)
    end
    nothing
end

function isdone(node::Node, x, options::Options)
    !node.main_node && return false

    options.done_callback(x, node.points) ||
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
