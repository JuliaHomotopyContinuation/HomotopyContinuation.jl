export monodromy_solve, realsolutions, nreal, GroupActions

################
# Group actions
################
"""
    GroupActions(actions::Function...)

Store a bunch of group actions `(f1, f2, f3, ...)`.
Each action has to return a tuple.
The actions are applied in the following sense
1) f1 is applied on the original solution `s`
2) f2 is applied on `s` and the results of 1
3) f3 is applied on `s` and the results of 1) and 2)
and so on

## Example
```julia-repl
julia> f1(s) = (s * s,);

julia> f2(s) = (2s, -s, 5s);

julia> f3(s) = (s + 1,);

julia> GroupActions(f1)(3)
(9,)

julia> GroupActions(f1,f2)(3)
(9, 18, -9, 45)

julia> GroupActions(f1,f2, f3)(3)
(9, 18, -9, 45, 10, 19, -8, 46)
```
"""
struct GroupActions{T<:Tuple}
    actions::T
end
GroupActions(::Nothing) = GroupActions(())
GroupActions(actions::GroupActions) = actions
GroupActions(actions::Function...) = GroupActions(actions)
GroupActions(actions) = GroupActions(actions...)

function (G::GroupActions)(solution::V) where V
    convert_if_necessary(apply_group_actions(G.actions, solution), solution)
end
convert_if_necessary(::Tuple{}, ::Any) = ()
convert_if_necessary(xs::NTuple{N, V}, ::V) where {N, V} = xs
convert_if_necessary(xs, ::V) where {V} = map(x -> convert(V, x), xs)

apply_group_actions(::Tuple{}, solution) = ()
# special case 1 function case
apply_group_actions(fs::Tuple{<:Function}, solution) = (fs[1])(solution)
# general case
function apply_group_actions(actions::Tuple, solution)
    f, rest = actions[1], Base.tail(actions)
    foldl(rest; init=f(solution)) do acc, f
        foldl(flatten, map(f, acc); init=acc)
    end
end
flatten(a::Tuple, b::Tuple) = (a..., b...)

#####################
# Monodromy Options #
#####################
const monodromy_options_allowed_keywords = [:identical_tol, :done_callback,
    :group_action,:group_actions, :group_action_on_all_nodes,
    :parameter_sampler, :equivalence_classes, :complex_conjugation, :target_solutions_count, :timeout,
    :minimal_number_of_solutions, :maximal_number_of_iterations_without_progress]

struct MonodromyOptions{F1<:Function, F2<:Tuple, F3<:Function}
    identical_tol::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
    group_action_on_all_nodes::Bool
    parameter_sampler::F3
    equivalence_classes::Bool
    complex_conjugation::Bool
    # stopping heuristic
    target_solutions_count::Int
    timeout::Float64
    minimal_number_of_solutions::Int
    maximal_number_of_iterations_without_progress::Int
end

function MonodromyOptions(isrealsystem;
    identical_tol::Float64=1e-6,
    done_callback=always_false,
    group_action=nothing,
    group_actions= group_action === nothing ? nothing : GroupActions(group_action),
    group_action_on_all_nodes=false,
    parameter_sampler=independent_normal,
    equivalence_classes=true,
    complex_conjugation=isrealsystem,
    # stopping heuristic
    target_solutions_count=nothing,
    timeout=float(typemax(Int)),
    minimal_number_of_solutions::Int=default_minimal_number_of_solutions(target_solutions_count),
    maximal_number_of_iterations_without_progress::Int=10)

    if group_actions isa GroupActions
       actions = group_actions
    else
       actions = GroupActions(group_actions)
    end

    MonodromyOptions(identical_tol, done_callback, actions,
        group_action_on_all_nodes, parameter_sampler, equivalence_classes, complex_conjugation,
        target_solutions_count == nothing ? typemax(Int) : target_solutions_count,
        float(timeout),
        minimal_number_of_solutions,
        maximal_number_of_iterations_without_progress)
end

default_minimal_number_of_solutions(::Nothing) = 2
function default_minimal_number_of_solutions(target_solutions_count::Int)
    div(target_solutions_count, 2)
end

always_false(x, sols) = false

has_group_actions(options::MonodromyOptions) = !(options.group_actions isa GroupActions{Tuple{}})


"""
    independent_normal(p::SVector{N, T}) where {N, T}

Sample a `SVector{N, T}` where each entries is drawn independently from the univariate normal distribution.
"""
function independent_normal(p::SVector{N, T}) where {N, T}
    @SVector randn(T, N)
end

##########################
## Monodromy Statistics ##
##########################

mutable struct MonodromyStatistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    nreal::Int
    nparametergenerations::Int
    nsolutions_development::Vector{Int}
end

MonodromyStatistics(nsolutions::Int) = MonodromyStatistics(0, 0, 0, 1, [nsolutions])
function MonodromyStatistics(solutions)
    stats = MonodromyStatistics(length(solutions))
    for s in solutions
        checkreal!(stats, s)
    end
    stats
end

Base.show(io::IO, S::MonodromyStatistics) = print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::MonodromyStatistics) = S

# update routines
function trackedpath!(stats::MonodromyStatistics, retcode)
    if retcode == PathTrackerStatus.success
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end

function checkreal!(stats::MonodromyStatistics, y)
    if isrealvector(y)
        stats.nreal +=1
    end
end

function generated_parameters!(stats::MonodromyStatistics, nsolutions::Int)
    stats.nparametergenerations += 1
    push!(stats.nsolutions_development, nsolutions)
end

function finished!(stats, nsolutions)
    push!(stats.nsolutions_development, nsolutions)
end

function n_loops_without_change(stats, nsolutions)
    k = 0
    for i in length(stats.nsolutions_development):-1:1
        if stats.nsolutions_development[i] != nsolutions
            return k
        end
        k += 1
    end
    return k
end

function n_solutions_current_loop(statistics, nsolutions)
    nsolutions - statistics.nsolutions_development[end]
end

#############################
# Loops and Data Structures #
#############################

export Triangle, Petal

#####################
# Nodes of the Loop
#####################

"""
    Node(p::SVector, x::AbstractVector; store_points=true, is_main_node=false)

Create a node with a parameter from the same type as `p` and expecting
points with tthe same type as `x`. If `stores_points` is `true` a `UniquePoints`
data structure is allocated to keep track of all known solutions of this node.
"""
struct Node{N, T, UP<:UniquePoints}
    p::SVector{N, T} # maybe just a Vector?
    points::Union{Nothing, UP}
    # Metadata / configuration
    main_node::Bool
end

function Node(p::SVector{N, T}, x::AbstractVector; store_points=true, is_main_node=false) where {N, T}
    uniquepoints = UniquePoints(typeof(x))
    if store_points == false
        Node{N, T, typeof(uniquepoints)}(p, nothing, is_main_node)
    else
        Node(p, uniquepoints, is_main_node)
    end
end
function Node(p::SVector{N, T}, node::Node{N,T,UP}; store_points=true, is_main_node=false) where {N, T, UP}
    if store_points == false
        Node{N, T, UP}(p, nothing, is_main_node)
    else
        Node{N, T, UP}(p, UniquePoints(UP), is_main_node)
    end
end

"""
    add!(node::Node, x; kwargs...)

Calls [`add!`](@ref) on the points of the Node.
"""
function add!(node::Node, x; kwargs...)
    if node.points === nothing
        false
    else
        add!(node.points, x; kwargs...)
    end
end

"""
    iscontained(node::Node, x; kwargs...)

Calls [`iscontained`](@ref) on the points of the Node.
"""
function iscontained(node::Node, x; kwargs...)
    if node.points === nothing
        false
    else
        iscontained(node.points, x; kwargs...)
    end
end


#####################
# Edges of the Loop
#####################
"""
    Edge(start::Int, target::Int; usegamma=false)

Create an edge between two nodes references by the index. `usegamma` refers
to an optional weight of two random complex numbers γ=(γ₁, γ₀). These are the same
gammas as in [`ParameterHomotopy`](@ref).

    Edge(edge::Edge)

Construct an indentical `Edge` as `edge`, but with newly samples γ (if applicable).
"""
struct Edge
    start::Int
    target::Int
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end
function Edge(start::Int, target::Int; usegamma=false)
    if usegamma
        γ = (randn(ComplexF64), randn(ComplexF64))
    else
        γ = nothing
    end
    Edge(start, target, γ)
end
function Edge(edge::Edge)
    if edge.γ === nothing
        edge
    else
        γ = (randn(ComplexF64), randn(ComplexF64))
        Edge(edge.start, edge.target, γ)
    end
end

#######################
# Loop data structure
#######################
"""
    Loop(p::SVector, x::AbstractVector{<:Number}, nnodes::Int, options::MonodromyOptions; usegamma=true)

Construct a loop using `nnodes` nodes with parameters of the type of `p` and solutions
of the type of `x`. `usegamma` refers to the use weights on the edges. See also
[`Edge`](@ref).

    Loop(p::SVector, x::AbstractVector, nnodes::Int, options::MonodromyOptions; usegamma=true)

Construct the loop and add all points in `x` to it.

    Loop(style::LoopStyle, p::SVector, x::AbstractVector, options::MonodromyOptions)

Construct a loop using the defined style `style` with parameters of the type of `p` and solutions
of the type of `x`.
"""
struct Loop{N<:Node}
    nodes::Vector{N}
    edges::Vector{Edge}
end

function Loop(p₁::SVector, x₁::AbstractVector{<:Number}, nnodes::Int, options::MonodromyOptions; usegamma=true)
    n₁ = Node(p₁, x₁, is_main_node = true)
    nodes = [n₁]
    store_points = options.group_action_on_all_nodes && has_group_actions(options)
    for i = 2:nnodes
        p = options.parameter_sampler(p₁)
        push!(nodes, Node(p, x₁, store_points=store_points))
    end

    loop = map(i -> Edge(i - 1, i; usegamma=usegamma), 2:nnodes)
    push!(loop, Edge(nnodes, 1; usegamma=usegamma))

    Loop(nodes, loop)
end
function Loop(p₁::SVector, x::AbstractVector, nnodes::Int, options::MonodromyOptions; kwargs...)
    loop = Loop(p₁, first(x), nnodes, options; kwargs...)
    for xᵢ ∈ x
        add!(loop.nodes[1], xᵢ; tol=options.identical_tol)
    end
    loop
end

"""
    solutions(loop::Loop)

Get the solutions of the loop.
"""
solutions(loop::Loop) = loop.nodes[1].points

"""
    nsolutions(loop::Loop)

Get the number solutions of the loop.
"""
nsolutions(loop::Loop) = length(solutions(loop))

"""
    mainnode(loop::Loop)

Get the main [`Node`](@ref) of the loop, i.e., the one with the original provided parameter.
"""
mainnode(loop::Loop) = loop.nodes[1]

nextedge(loop::Loop, edge::Edge) = loop.edges[edge.target]

##############
# Loop styles
##############
"""
    LoopStyle

Abstract type defining a style of a loop.
"""
abstract type LoopStyle end


"""
    Triangle(;weights=true)

A triangle is a loop consisting of the main node and two addtional nodes.
If `weights` is true the edges are equipped with additional random weights.
Note that this is usually only necessary for real parameters.
"""
struct Triangle <: LoopStyle
    useweights::Bool
end
Triangle(;useweights=true) = Triangle(useweights)

function Loop(strategy::Triangle, p::SVector, x::AbstractVector, options::MonodromyOptions)
    Loop(p, x, 3, options, usegamma=strategy.useweights)
end

"""
    Petal()

A petal is a loop consisting of the main node and one other node connected
by two edges with different random weights.
"""
struct Petal <: LoopStyle end
function Loop(strategy::Petal, p::SVector, x::AbstractVector, options::MonodromyOptions)
    Loop(p, x, 2, options, usegamma=true)
end


"""
    regenerate!(loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics)

Regenerate all random parameters in the loop in order to introduce a new monodromy action.
"""
function regenerate!(loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics)
    main = mainnode(loop)

    # The first node is the main node and doesn't get touched
    for i ∈ 2:length(loop.nodes)
        loop.nodes[i] = Node(options.parameter_sampler(main.p), loop.nodes[i])
    end
    loop.edges .= Edge.(loop.edges)
    generated_parameters!(stats, length(main.points)) # bookkeeping
end


"""
    track(tracker, x::AbstractVector, edge::Edge, loop::Loop, stats::MonodromyStatistics)

Track `x` along the edge `edge` in the loop `loop` using `tracker`. Record statistics
in `stats`.
"""
function track(tracker::PathTracker, x::AbstractVector, edge::Edge, loop::Loop, stats::MonodromyStatistics)
    set_parameters!(tracker, edge, loop)
    track(tracker, x, stats)
end
function track(tracker::PathTracker, x::AbstractVector, stats::MonodromyStatistics)
    retcode = track!(tracker, x, 1.0, 0.0)
    trackedpath!(stats, retcode)
    retcode
end

"""x
    set_parameters!(tracker::PathTracker, e::Edge, loop::Loop)

Setup the parameters in the ParameterHomotopy in `tracker` to fit the edge `e`.
"""
function set_parameters!(tracker::PathTracker, e::Edge, loop::Loop)
    H = basehomotopy(tracker.homotopy)
    if !(H isa ParameterHomotopy)
        error("Base homotopy is not a ParameterHomotopy")
    end
    p = (loop.nodes[e.start].p, loop.nodes[e.target].p)
    set_parameters!(H, p, e.γ)
end


#############
## Results ##
#############
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
    mapresults(f, result::MonodromyResult; onlyreal=false, realtol=1e-6)

Apply the function `f` to all entries of `MonodromyResult` for which the given conditions apply.

## Example
```julia
# This gives us all solutions considered real (but still as a complex vector).
realsolutions = mapresults(solution, R, onlyreal=true)
```
"""
function mapresults(f, R::MonodromyResult;
    onlyreal=false, realtol=1e-6)
    [f(r) for r in R.solutions if
        (!onlyreal || isrealvector(r, realtol))]
end

"""
    solutions(result::MonodromyResult; onlyreal=false, realtol=1e-6)

Return all solutions (as `SVector`s) for which the given conditions apply.

## Example
```julia
realsolutions = solutions(R, onlyreal=true)
```
"""
function solutions(R::MonodromyResult; kwargs...)
    mapresults(identity, R; kwargs...)
end

"""
    nsolutions(result::MonodromyResult)

Returns the number solutions of the `result`.
"""
nsolutions(res::MonodromyResult) = length(res.solutions)

"""
    realsolutions(res::MonodromyResult; tol=1e-6)

Returns the solutions of `res` whose imaginary part has norm less than 1e-6.
"""
function realsolutions(res::MonodromyResult; tol=1e-6)
    map(r -> real.(r), filter(r -> LinearAlgebra.norm(imag.(r)) < tol, res.solutions))
end

"""
    nreal(res::MonodromyResult; tol=1e-6)

Counts how many solutions of `res` have imaginary part norm less than 1e-6.
"""
function nreal(res::MonodromyResult; tol=1e-6)
    count(r -> LinearAlgebra.norm(imag.(r)) < tol, res.solutions)
end


#####################
## monodromy solve ##
#####################
"""
MonodromyCache{FT<:FixedHomotopy, Tracker<:PathTracker, NC<:NewtonCache}

Cache for monodromy loops.
"""
struct MonodromyCache{FT<:FixedHomotopy, Tracker<:PathTracker, NC<:NewtonCache, AV<:AbstractVector}
    F::FT
    tracker::Tracker
    newton_cache::NC
    out::AV
end


"""
    monodromy_solve(F, sols, p; parameters=..., options..., pathtrackerkwargs...)

Solve a polynomial system `F(x;p)` with specified parameters and initial solutions `sols`
by monodromy techniques. This makes loops in the parameter space of `F` to find new solutions.

## Options
* `target_solutions_count=nothing`: The computations are stopped if this number of solutions is reached.
* `done_callback=always_false`: A callback to end the computation early. This function takes 2 arguments. The first one is the new solution `x` and the second one are all current solutions (including `x`). Return `true` if the compuation is done.
* `maximal_number_of_iterations_without_progress::Int=10`: The maximal number of iterations (i.e. loops generated) without any progress.
* `group_action=nothing`: A function taking one solution and returning other solutions if there is a constructive way to obtain them, e.g. by symmetry.
* `strategy`: The strategy used to create loops. If `F` only depends linearly on `p` this will be [`Petal`](@ref). Otherwise this will be [`Triangle`](@ref) with weights if `F` is a real system.
* `showprogress=true`: Enable a progress meter.
* `identical_tol::Float64=1e-6`: The tolerance with which it is decided whether two solutions are identical.
* `group_actions=nothing`: If there is more than one group action you can use this to chain the application of them. For example if you have two group actions `foo` and `bar` you can set `group_actions=[foo, bar]`. See [`GroupActions`](@ref) for details regarding the application rules.
* `group_action_on_all_nodes=false`: By default the group_action(s) are only applied on the solutions with the main parameter `p`. If this is enabled then it is applied for every parameter `q`.
* `parameter_sampler=independent_normal`: A function taking the parameter `p` and returning a new random parameter `q`. By default each entry of the parameter vector is drawn independently from the univariate normal distribution.
* `equivalence_classes=true`: This only applies if there is at least one group action supplied. We then consider two solutions in the same equivalence class if we can transform one to the other by the supplied group actions. We only track one solution per equivalence class.
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
        parameters=throw(ArgumentError("You need to provide `parameters=...` to monodromy")),
        strategy=default_strategy(F, parameters, p₀),
        showprogress=true,
        kwargs...) where {TC, TP, NParams, NVars}

    if length(p₀) ≠ length(parameters)
        throw(ArgumentError("Number of provided parameters doesn't match the length of initially provided parameter `p₀`."))
    end

    p₀ = convert(SVector{NParams, promote_type(Float64, TP)}, p₀)

    optionskwargs, restkwargs = splitkwargs(kwargs, monodromy_options_allowed_keywords)
    options = begin
        isrealsystem = TC <: Real && TP <: Real
        MonodromyOptions(isrealsystem; optionskwargs...)
    end

    #assemble
    loop = Loop(strategy, p₀, startsolutions, options)

    tracker = pathtracker(F, startsolutions;
                          parameters=parameters, p₁=p₀, p₀=p₀, restkwargs...)
    statistics = MonodromyStatistics(solutions(loop))

    if affine_tracking(tracker)
        HC = HomotopyWithCache(tracker.homotopy, tracker.state.x, 1.0)
        F₀ = FixedHomotopy(HC, 0.0)
    else
        # Force affine newton method
        patch_state = state(EmbeddingPatch(), tracker.state.x)
        HC = HomotopyWithCache(PatchedHomotopy(tracker.homotopy, patch_state), tracker.state.x, 1.0)
        F₀ = FixedHomotopy(HC, 0.0)
    end
    # construct cache
    newton_cache = NewtonCache(F₀, tracker.state.x)
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

function default_strategy(F::Vector{<:MP.AbstractPolynomialLike{TC}}, parameters, p₀::AbstractVector{TP}) where {TC,TP}
    # If F depends only linearly on the parameters a petal is sufficient
    if all(f -> last(minmaxdegree(f, parameters)) ≤ 1, F)
        Petal()
    # For a real system we should introduce some weights to avoid the discriminant
    elseif TP <: Real && TC <: Real
        Triangle(useweights=true)
    else
        Triangle(useweights=false)
    end
end

# convert vector of vectors to vector of svectors
static_solutions(V::Vector) = static_solutions(V, Val(length(V[1])))
function static_solutions(V::Vector, ::Val{N}) where {N}
    map(v -> complex.(float.(SVector{N}(v))), V)
end
function static_solutions(V::Vector{<:AbstractVector{<:Complex{<:AbstractFloat}}}, ::Val{N}) where {N}
    SVector{N}.(V)
end


#################
## Actual work ##
#################

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
        elseif retcode == :invalid_startvalue
            return :invalid_startvalue
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
        status = process!(queue, job, C, loop, options, stats, progress)
        if status == :done
            return :done
        elseif status == :invalid_startvalue
            return :invalid_startvalue
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
    # We distinguish solutions which have a distance larger than identical_tol
    # But due to the numerical error in the evaluation of the distance, we need to be a little bit
    # carfule. Therefore, we require that the solutions should be one magnitude closer to
    # the true solutions as necessary
    tol = 0.1 * options.identical_tol
    result = newton!(C.out, C.F, ŷ, euclidean_norm, C.newton_cache,
                tol=tol, miniters=1, maxiters=3, simplified_last_step=false)

    if result.retcode == converged
        return affine_chart(x, C.out)
    else
        return nothing
    end
end

affine_chart(x::SVector, y::PVector) = ProjectiveVectors.affine_chart!(x, y)
affine_chart(x::SVector{N, T}, y::AbstractVector) where {N, T} = SVector{N,T}(y)

function process!(queue::Vector{<:Job}, job::Job, C::MonodromyCache, loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics, progress)
    retcode = track(C.tracker, job.x, job.edge, loop, stats)
    if retcode ≠ PathTrackerStatus.success
        if retcode == PathTrackerStatus.terminated_invalid_startvalue && stats.ntrackedpaths == 0
            return :invalid_startvalue
        end
        return :incomplete
    end

    node = loop.nodes[job.edge.target]

    if node.main_node
        y = verified_affine_vector(C, currx(C.tracker), job.x, options)
        # is the solution at infinity?
        if y === nothing
            return :incomplete
        end
    else
        y = affine_chart(job.x, currx(C.tracker))
    end


    if !iscontained(node, y; tol=options.identical_tol)
        if options.equivalence_classes
            if equivalence_class_contained(node, y, options)
                return :incomplete
            end
        end

        add!(node, y; tol=options.identical_tol)
        # If we are on the main node check whether we have a real root.
        node.main_node && checkreal!(stats, y)
        # Check if we are done
        node.main_node && isdone(node, y, options) && return :done

        next_edge = nextedge(loop, job.edge)
        push!(queue, Job(y, next_edge))

        # Check for complex conjugate solution
        if options.complex_conjugation && node.main_node
            ȳ = conj.(y)
            if !equivalence_class_contained(node, ȳ, options)
                add!(node, ȳ; tol=options.identical_tol)
                # Check if we are done
                isdone(node, ȳ, options) && return :done
                # Schedule new job
                push!(queue, Job(ȳ, next_edge))
            end
        end

        # Handle group actions
        # 1) We only need to go through group actions if we do not compute in equivalence classes
        # 2) Things are setup up such that for nodes where we want to apply
        #    group actions `node.points !== nothing`
        if !options.equivalence_classes && node.points !== nothing
            for yᵢ in options.group_actions(y)
                if add!(node, yᵢ; tol=options.identical_tol)
                    node.main_node && checkreal!(stats, yᵢ)
                    # Check if we are done
                    node.main_node && isdone(node, yᵢ, options) && return :done

                    push!(queue, Job(yᵢ, next_edge))
                end
            end
        end
    end
    return :incomplete
end

function equivalence_class_contained(node, y, options)
    for yᵢ in options.group_actions(y)
        if iscontained(node, yᵢ, tol=options.identical_tol)
            # equivalence class already existing
            return true
        end
    end
    false
end

function update_progress!(::Nothing, loop::Loop, statistics::MonodromyStatistics; finish=false)
    nothing
end
function update_progress!(progress, loop::Loop, statistics::MonodromyStatistics; finish=false)
    nsolutions = length(solutions(loop))
    ProgressMeter.update!(progress, nsolutions, showvalues=(
        ("# paths tracked", statistics.ntrackedpaths),
        ("# loops generated", statistics.nparametergenerations),
        ("# loops without change", n_loops_without_change(statistics, nsolutions)),
        ("# solutions in current loop", n_solutions_current_loop(statistics, nsolutions)),
        ("# real solutions", statistics.nreal),
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
