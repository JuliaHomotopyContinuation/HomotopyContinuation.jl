export monodromy_solve, realsolutions, nreal, GroupActions, complex_conjugation

include("monodromy/group_actions.jl")
include("monodromy/options.jl")
include("monodromy/statistics.jl")
include("monodromy/loop.jl")

struct MonodromyResult{N, T}
    returncode::Symbol
    solutions::Vector{SVector{N, T}}
    statistics::MonodromyStatistics
end

Base.iterate(R::MonodromyResult) = iterate(R.solutions)
Base.iterate(R::MonodromyResult, state) = iterate(R.solutions, state)

Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MonodromyResult) = x
function Base.show(io::IO, result::MonodromyResult{N, T}) where {N, T}
    println(io, "MonodromyResult")
    println(io, "==================================")
    println(io, "• $(nsolutions(result)) solutions ($(nreal(result)) real)")
    println(io, "• return code → $(result.returncode)")
    println(io, "• $(result.statistics.ntrackedpaths) tracked paths")
end


TreeViews.hastreeview(::MonodromyResult) = true
TreeViews.numberofnodes(::MonodromyResult) = 4
TreeViews.treelabel(io::IO, x::MonodromyResult, ::MIME"application/prs.juno.inline") =
    print(io, "<span class=\"syntax--support syntax--type syntax--julia\">MonodromyResult</span>")

function TreeViews.nodelabel(io::IO, x::MonodromyResult, i::Int, ::MIME"application/prs.juno.inline")
    if i == 1
        print(io, "Solutions")
    elseif i == 2
            print(io, "Real solutions")
    elseif i == 3
        print(io, "Return Code")
    elseif i == 4
        print(io, "Tracked Paths")
    end
end
function TreeViews.treenode(r::MonodromyResult, i::Integer)
    if i == 1
        return r.solutions
    elseif i == 2
        return realsolutions(r)
    elseif i == 3
        return r.returncode
    elseif i == 4
        return r.statistics.ntrackedpaths
    end
    missing
end


"""
    solutions(result::MonodromyResult)

Returns the solutions of the `result`.
"""
solutions(res::MonodromyResult) = res.solutions

"""
    nsolutions(result::MonodromyResult)

Returns the number solutions of the `result`.
"""
nsolutions(res::MonodromyResult) = length(res.solutions)

"""
    realsolutions(res::MonodromyResult; tol = 1e-6)

Returns the solutions of `res` whose imaginary part has norm < 1e-6.
"""
function realsolutions(res::MonodromyResult; tol=1e-6)
    map(r -> real.(r), filter(r -> LinearAlgebra.norm(imag.(r)) < tol, res.solutions))
end

"""
    nreal(res::MonodromyResult; tol = 1e-6)

Counts how many solutions of `res` have imaginary part norm < 1e-6.
"""
function nreal(res::MonodromyResult; tol=1e-6)
    count(r -> LinearAlgebra.norm(imag.(r)) < tol, res.solutions)
end

"""
MonodromyCache{FT<:FixedHomotopy, Tracker<:PathTracker, NC<:NewtonCache}

Cache for monodromy loops.
"""
struct MonodromyCache{FT<:FixedHomotopy, Tracker<:PathTracker, NC<:NewtonCache, T<:Number, N}
    F::FT
    tracker::Tracker
    newton_cache::NC
    out::ProjectiveVectors.PVector{T,N}
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
* `strategy`: The strategy used to create loops. By default this will be `Triangle` with weights if `F` is a real system.
* `showprogress=true`: Enable a progress meter.
* `accuracy::Float64=1e-6`: The tolerance with which it is decided whether two solutions are identical.
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
        scale_system=true,
        showprogress=true,
        kwargs...) where {TC, TP, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        error("Number of provided parameters doesn't match the length of initially provided parameter `p₀`.")
    end

    p₀ = convert(SVector{NParams, promote_type(Float64, TP)}, p₀)

    optionskwargs, restkwargs = splitkwargs(kwargs, options_allowed_keywords)
    options = begin
        isrealsystem = TC <: Real && TP <: Real
        MonodromyOptions(isrealsystem; optionskwargs...)
    end

    #assemble
    loop = Loop(strategy, p₀, startsolutions, options)
    if scale_system
        f = F ./ coefficient_norm.(map(fi -> MP.subs(fi, parameters=>p₀), F))
    else
        f = F
    end

    tracker = pathtracker(
        f, startsolutions; parameters=parameters, p₁=p₀, p₀=p₀, restkwargs...)
    statistics = MonodromyStatistics(solutions(loop))

    # affine newton methods
    patch_state = state(EmbeddingPatch(), tracker.state.x)
    HC = HomotopyWithCache(PatchedHomotopy(tracker.homotopy, patch_state), tracker.state.x, 1.0)
    F₀ = FixedHomotopy(HC, 0.0)
    newton_cache = NewtonCache(F₀, tracker.state.x)

    # construct cache
    C =  MonodromyCache(F₀, tracker, newton_cache, copy(tracker.state.x))


    # solve
    retcode = :not_assigned
    if showprogress
        progress = ProgressMeter.ProgressUnknown("Solutions found:")
    else
        progress = nothing
    end
    try
        retcode = monodromy_solve!(loop, C, options, statistics, progress)
    catch e
        if (e isa InterruptException)
            retcode = :interrupt
        else
            rethrow(e)
        end
    end
    finished!(statistics, nsolutions(loop))
    MonodromyResult(retcode, points(solutions(loop)), statistics)
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


function monodromy_solve!(loop::Loop, C::MonodromyCache, options::MonodromyOptions,
    stats::MonodromyStatistics, progress)

    t₀ = time_ns()
    iterations_without_progress = 0 # stopping heuristic
    # intialize job queue
    queue = map(x -> Job(x, loop.edges[1]), solutions(loop))

    n = nsolutions(loop)
    while n < options.target_solutions_count
        retcode = empty_queue!(queue, loop, C, options, t₀, stats, progress)

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

function empty_queue!(queue, loop::Loop, C::MonodromyCache, options::MonodromyOptions,
        t₀::UInt, stats::MonodromyStatistics, progress)
    while !isempty(queue)
        job = pop!(queue)
        if process!(queue, job, C, loop, options, stats, progress) == :done
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

function verified_affine_vector(C::MonodromyCache, ŷ, x, options)
    result = newton!(C.out, C.F, ŷ, options.accuracy, 3, true, 1.0, C.newton_cache)

    if result.retcode == converged
        return ProjectiveVectors.affine_chart!(x, C.out)
    else
        return nothing
    end
end

function process!(queue::Vector{<:Job}, job::Job, C::MonodromyCache, loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics, progress)
    retcode = track(C.tracker, job.x, job.edge, loop, stats)
    if retcode ≠ PathTrackerStatus.success
        return :incomplete
    end

    if job.edge.target == 1
        y = verified_affine_vector(C, currx(C.tracker), job.x, options)
        #is the solution at infinity?
        if y === nothing
            return :incomplete
        end
    else
        y = ProjectiveVectors.affine_chart!(job.x, currx(C.tracker))
    end

    if job.edge.target == 1
        #is the solution real?
        checkreal!(stats, y)
    end

    node = loop.nodes[job.edge.target]
    if !iscontained(node, y, tol=options.accuracy)
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
                if !iscontained(node, yᵢ, tol=options.accuracy)
                    unsafe_add!(node, yᵢ)
                    if job.edge.target == 1
                        checkreal!(stats, y)
                    end
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

function update_progress!(::Nothing, loop::Loop, statistics::MonodromyStatistics; finish=false)
    nothing
end
function update_progress!(progress, loop::Loop, statistics::MonodromyStatistics; finish=false)
    ProgressMeter.update!(progress, length(solutions(loop)), showvalues=(
        ("# paths tracked ", statistics.ntrackedpaths),
        ("# loops generated ", statistics.nparametergenerations),
        ("# real solutions ", statistics.nreal)
    ))
    if finish
        ProgressMeter.finish!(progress)
    end
    nothing
end

function isdone(node::Node, x, options::MonodromyOptions)
    !node.main_node && return false

    options.done_callback(x, node.points) ||
    length(node.points) ≥ options.target_solutions_count
end

function regenerate_loop_and_schedule_jobs!(queue, loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics)
    sols = solutions(loop)
    # create a new loop by regenerating the parameters (but don't touch our
    # main node)
    regenerate!(loop, options, stats)
    for x in sols
        push!(queue, Job(x, loop.edges[1]))
    end
    nothing
end
