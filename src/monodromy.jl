export monodromy_solve,
       MonodromyResult,
       real_solutions,
       is_success,
       is_heuristic_stop,
       nreal,
       parameters,
       verify_solution_completeness,
       solution_completeness_witnesses


#####################
# Monodromy Options #
#####################
const monodromy_options_supported_keywords = [
    :distance,
    :identical_tol,
    :done_callback,
    :group_action,
    :group_actions,
    :group_action_on_all_nodes,
    :parameter_sampler,
    :equivalence_classes,
    :complex_conjugation,
    :check_startsolutions,
    :target_solutions_count,
    :timeout,
    :min_solutions,
    :max_loops_no_progress,
    :reuse_loops,
]

struct MonodromyOptions{F<:Function,F1<:Function,F2<:Tuple,F3<:Function}
    distance_function::F
    identical_tol::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
    group_action_on_all_nodes::Bool
    parameter_sampler::F3
    equivalence_classes::Bool
    complex_conjugation::Bool
    check_startsolutions::Bool
    # stopping heuristic
    target_solutions_count::Int
    timeout::Float64
    min_solutions::Int
    max_loops_no_progress::Int
    reuse_loops::Symbol
end

function MonodromyOptions(
    is_real_system::Bool,
    accuracy::Float64,
    distance;
    identical_tol::Float64 = sqrt(accuracy),
    done_callback = always_false,
    group_action = nothing,
    group_actions = group_action === nothing ? nothing : GroupActions(group_action),
    group_action_on_all_nodes = false,
    parameter_sampler = independent_normal,
    equivalence_classes = true,
    complex_conjugation = is_real_system,
    check_startsolutions = true,
    # stopping heuristic
    target_solutions_count = nothing,
    timeout = float(typemax(Int)),
    min_solutions::Int = default_min_solutions(target_solutions_count),
    max_loops_no_progress::Int = 10,
    reuse_loops::Symbol = :all,
)

    if group_actions isa GroupActions
        actions = group_actions
    else
        if group_actions === nothing
            equivalence_classes = false
        end
        actions = GroupActions(group_actions)
    end


    MonodromyOptions(
        distance,
        identical_tol,
        done_callback,
        actions,
        group_action_on_all_nodes,
        parameter_sampler,
        equivalence_classes,
        complex_conjugation,
        check_startsolutions,
        target_solutions_count === nothing ? typemax(Int) : target_solutions_count,
        float(timeout),
        min_solutions,
        max_loops_no_progress,
        reuse_loops,
    )
end

default_min_solutions(::Nothing) = 1
function default_min_solutions(target_solutions_count::Int)
    div(target_solutions_count, 2)
end

always_false(x, sols) = false

has_group_actions(options::MonodromyOptions) =
    !(options.group_actions isa GroupActions{Tuple{}})


"""
    independent_normal(p::AbstractVector{T}) where {T}

Sample a vector where each entries is drawn independently from the univariate normal distribution.
"""
independent_normal(p::SVector{N,T}) where {N,T} = @SVector randn(T, N)
independent_normal(p::AbstractVector{T}) where {T} = randn(T, length(p))

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
        if is_real_vector(s)
            stats.nreal += 1
        end
    end
    stats
end

Base.show(io::IO, S::MonodromyStatistics) = print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::MonodromyStatistics) = S

# update routines
function trackedpath!(stats::MonodromyStatistics, retcode)
    if is_success(retcode)
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end

function generated_parameters!(stats::MonodromyStatistics, nsolutions::Int)
    stats.nparametergenerations += 1
    push!(stats.nsolutions_development, nsolutions)
end

function finished!(stats, nsolutions)
    push!(stats.nsolutions_development, nsolutions)
end

function n_completed_loops_without_change(stats, nsolutions)
    k = 0
    for i = length(stats.nsolutions_development):-1:1
        if stats.nsolutions_development[i] != nsolutions
            return max(k - 1, 0)
        end
        k += 1
    end
    return max(k - 1, 0)
end

function n_solutions_current_loop(statistics, nsolutions)
    nsolutions - statistics.nsolutions_development[end]
end

#############################
# Loops and Data Structures #
#############################

export Triangle, Petal

#######################
# Loop data structure
#######################

struct LoopEdge
    p₁::Vector{ComplexF64}
    p₀::Vector{ComplexF64}
    weights::Union{Nothing,NTuple{2,ComplexF64}}
end

function LoopEdge(p₁, p₀; weights = false)
    if weights
        γ = (randn(ComplexF64), randn(ComplexF64))
    else
        γ = nothing
    end
    LoopEdge(convert(Vector{ComplexF64}, p₁), convert(Vector{ComplexF64}, p₀), γ)
end

struct MonodromyLoop
    edges::Vector{LoopEdge}
end

function MonodromyLoop(base_p, nnodes::Int, options::MonodromyOptions; weights = true)
    edges = LoopEdge[]
    p₁ = base_p
    for i = 2:nnodes
        p₀ = options.parameter_sampler(p₁)
        push!(edges, LoopEdge(p₁, p₀; weights = weights))
        p₁ = p₀
    end
    push!(edges, LoopEdge(p₁, base_p; weights = weights))
    MonodromyLoop(edges)
end

##############
# Loop styles
##############
"""
    LoopStyle

Abstract type defining a style of a loop.
"""
abstract type LoopStyle end


"""
    Triangle(;useweights=true)

A triangle is a loop consisting of the main node and two addtional nodes.
If `weights` is true the edges are equipped with additional random weights.
Note that this is usually only necessary for real parameters.
"""
struct Triangle <: LoopStyle
    useweights::Bool
end
Triangle(; useweights = true) = Triangle(useweights)

function MonodromyLoop(strategy::Triangle, p, options::MonodromyOptions)
    MonodromyLoop(p, 3, options, weights = strategy.useweights)
end

"""
    Petal()

A petal is a loop consisting of the main node and one other node connected
by two edges with different random weights.
"""
struct Petal <: LoopStyle end
function MonodromyLoop(strategy::Petal, p, options::MonodromyOptions)
    MonodromyLoop(p, 2, options, weights = true)
end


"""
    regenerate!(loop::MonodromyLoop, options::MonodromyOptions, stats::MonodromyStatistics)

Regenerate all random parameters in the loop in order to introduce a new monodromy action.
"""
function regenerate!(
    loop::MonodromyLoop,
    options::MonodromyOptions,
    stats::MonodromyStatistics,
)
    main = mainnode(loop)

    # The first node is the main node and doesn't get touched
    store_points = options.group_action_on_all_nodes && has_group_actions(options)
    for i ∈ 2:length(loop.nodes)
        loop.nodes[i] = Node(
            options.parameter_sampler(main.p),
            loop.nodes[i],
            options;
            store_points = store_points,
            is_main_node = false,
        )
    end
    loop.edges .= Edge.(loop.edges)
    generated_parameters!(stats, length(main.points)) # bookkeeping
end


"""
    track(tracker, x, edge::LoopEdge, loop::MonodromyLoop, stats::MonodromyStatistics)

Track `x` along the edge `edge` in the loop `loop` using `tracker`. Record statistics
in `stats`.
"""
function track(
    tracker::CoreTracker,
    x::AbstractVector,
    loop::MonodromyLoop,
    ntrackedpaths::Base.RefValue{Int},
    ntrackingfailures::Base.RefValue{Int},
)
    H = basehomotopy(tracker.homotopy)
    local retcode::CoreTrackerStatus.states
    y = x
    for e in loop.edges
        set_parameters!(H, (e.p₁, e.p₀), e.weights)
        retcode = track!(tracker, y)
        if is_success(retcode)
            ntrackedpaths[] += 1
        else
            ntrackingfailures[] += 1
        end
        is_success(retcode) || break
        y = current_x(tracker)
    end
    retcode
end


#############
## Results ##
#############
"""
    MonodromyResult

The monodromy result contains the result of the `monodromy_solve` computation.
"""
struct MonodromyResult{AV<:AbstractVector,T}
    returncode::Symbol
    solutions::Vector{AV}
    parameters::Vector{T}
    statistics::MonodromyStatistics
    equivalence_classes::Bool
    seed::Int
end

Base.iterate(R::MonodromyResult) = iterate(R.solutions)
Base.iterate(R::MonodromyResult, state) = iterate(R.solutions, state)

Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MonodromyResult) = x
function Base.show(io::IO, result::MonodromyResult{N,T}) where {N,T}
    println(io, "MonodromyResult")
    println(io, "==================================")
    if result.equivalence_classes
        println(
            io,
            "• ",
            nsolutions(result),
            " classes of solutions (modulo group action) (",
            nreal(result),
            ") real)",
        )
    else
        println(io, "• $(nsolutions(result)) solutions ($(nreal(result)) real)")
    end
    println(io, "• return code → $(result.returncode)")
    println(io, "• $(result.statistics.ntrackedpaths) tracked paths")
    println(io, "• seed → $(result.seed)")
end


TreeViews.hastreeview(::MonodromyResult) = true
TreeViews.numberofnodes(::MonodromyResult) = 6
TreeViews.treelabel(io::IO, x::MonodromyResult, ::MIME"application/prs.juno.inline") =
    print(
        io,
        "<span class=\"syntax--support syntax--type syntax--julia\">MonodromyResult</span>",
    )

function TreeViews.nodelabel(
    io::IO,
    x::MonodromyResult,
    i::Int,
    ::MIME"application/prs.juno.inline",
)
    if i == 1
        if x.equivalence_classes
            print(io, "$(nsolutions(x)) classes of solutions (modulo group action)")
        else
            print(io, "$(nsolutions(x)) solutions")
        end
    elseif i == 2
        if x.equivalence_classes
            print(io, "$(nreal(x)) classes of real solutions")
        else
            print(io, "$(nreal(x)) real solutions")
        end
    elseif i == 3
        print(io, "Return code")
    elseif i == 4
        print(io, "Statistics")
    elseif i == 5
        print(io, "Parameters")
    elseif i == 6
        print(io, "Seed")
    end
end
function TreeViews.treenode(r::MonodromyResult, i::Integer)
    if i == 1
        return r.solutions
    elseif i == 2
        return real_solutions(r)
    elseif i == 3
        return r.returncode
    elseif i == 4
        return r.statistics
    elseif i == 5
        return r.parameters
    elseif i == 6
        return r.seed
    end
    missing
end

"""
    is_success(result::MonodromyResult)

Returns true if the monodromy computation achieved its target solution count.
"""
is_success(result::MonodromyResult) = result.returncode == :success

"""
    is_heuristic_stop(result::MonodromyResult)

Returns true if the monodromy computation stopped due to the heuristic.
"""
is_heuristic_stop(result::MonodromyResult) = result.returncode == :heuristic_stop

"""
    mapresults(f, result::MonodromyResult; only_real=false, real_tol=1e-6)

Apply the function `f` to all entries of `MonodromyResult` for which the given conditions apply.

## Example
```julia
# This gives us all solutions considered real (but still as a complex vector).
real_solutions = mapresults(solution, R, only_real=true)
```
"""
function mapresults(f, R::MonodromyResult; only_real = false, real_tol = 1e-6)
    [f(r) for r in R.solutions if (!only_real || is_real_vector(r, real_tol))]
end

"""
    solutions(result::MonodromyResult; only_real=false, real_tol=1e-6)

Return all solutions (as `SVector`s) for which the given conditions apply.

## Example
```julia
real_solutions = solutions(R, only_real=true)
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
    real_solutions(res::MonodromyResult; tol=1e-6)

Returns the solutions of `res` whose imaginary part has norm less than 1e-6.
"""
function real_solutions(res::MonodromyResult; tol = 1e-6)
    map(r -> real_vector(r), filter(r -> is_real_vector(r, tol), res.solutions))
end

"""
    nreal(res::MonodromyResult; tol=1e-6)

Counts how many solutions of `res` have imaginary part norm less than 1e-6.
"""
function nreal(res::MonodromyResult; tol = 1e-6)
    count(r -> is_real_vector(r, tol), res.solutions)
end

"""
    parameters(r::MonodromyResult)

Return the parameters corresponding to the given result `r`.
"""
parameters(r::MonodromyResult) = r.parameters

#####################
## monodromy solve ##
#####################
struct MonodromySolver{
    T<:Number,
    UP<:UniquePoints,
    MO<:MonodromyOptions,
    Tracker<:CoreTracker,
    Prob<:Problem,
}
    parameters::Vector{T}
    solutions::UP
    solutions_lock::ReentrantLock
    loops::Vector{MonodromyLoop}
    options::MO
    statistics::MonodromyStatistics
    trackers::Vector{Tracker}
    problem::Prob
end

function MonodromySolver(
    F::Inputs,
    solution::AbstractVector{<:Number},
    p₀::AbstractVector{TP};
    kwargs...,
) where {TP}
    MonodromySolver(F, [solution], p₀; kwargs...)
end
function MonodromySolver(
    F::Inputs,
    startsolutions::Vector{<:AbstractVector{<:Number}},
    p::AbstractVector{TP};
    distance = default_distance(eltype(startsolutions)),
    parameters = nothing,
    strategy = nothing,
    show_progress = true,
    showprogress = nothing, #deprecated
    accuracy = 1e-6, # set a higher refinement_accuracy than the default
    precision = PRECISION_ADAPTIVE,
    kwargs...,
) where {TP,NVars}

    @deprecatekwarg showprogress show_progress

    if parameters !== nothing && length(p) ≠ length(parameters)
        throw(ArgumentError("Number of provided parameters doesn't match the length of initially provided parameter `p₀`."))
    end

    x₀ = complex.(first(startsolutions))
    p₀ = Vector{promote_type(Float64, TP)}(p)

    optionskwargs, restkwargs = splitkwargs(kwargs, monodromy_options_supported_keywords)

    supported, rest = splitkwargs(restkwargs, problem_startsolutions_supported_keywords)
    prob, _ = problem_startsolutions(
        F,
        [x₀];
        parameters = parameters,
        generic_parameters = p₀,
        supported...,
    )
    tracker = CoreTracker(
        prob,
        x₀;
        accuracy = accuracy,
        precision = PRECISION_ADAPTIVE,
        rest...,
    )
    # Check whether homotopy is real
    is_real_system = numerically_check_real(tracker.homotopy, x₀)
    options = MonodromyOptions(is_real_system, accuracy, distance; optionskwargs...)
    # construct UniquePoints

    uniquepoints = UniquePoints(
        typeof(x₀),
        distance;
        group_actions = options.equivalence_classes ? options.group_actions : nothing,
        check_real = true,
    )

    # add only unique points that are true solutions
    for s in startsolutions
        if !options.check_startsolutions || is_valid_start_value(tracker, s, 0.0)
            add!(uniquepoints, s; tol = options.identical_tol)
        end
    end

    statistics = MonodromyStatistics(uniquepoints)

    if strategy === nothing
        strategy = default_strategy(F, parameters, p; is_real_system = is_real_system)
    end

    # construct Loop
    loops = [MonodromyLoop(strategy, p₀, options)]

    MonodromySolver(
        p₀,
        uniquepoints,
        ReentrantLock(),
        loops,
        options,
        statistics,
        [tracker],
        prob,
    )
end

"""
    solutions(MS::MonodromySolver)

Get the solutions of the loop.
"""
solutions(MS::MonodromySolver) = MS.solutions.points

"""
    nsolutions(loop::MonodromyLoop)

Get the number solutions of the loop.
"""
nsolutions(MS::MonodromySolver) = length(solutions(MS))


"""
    monodromy_solve(F, [sols, p]; parameters=..., options..., pathtrackerkwargs...)

Solve a polynomial system `F(x;p)` with specified parameters and initial solutions `sols`
by monodromy techniques. This makes loops in the parameter space of `F` to find new solutions.
If `F` the parameters `p` only occur *linearly* in `F` it is eventually possible to compute
a *start pair* ``(x₀, p₀)`` automatically. In this case `sols` and `p` can be omitted and
the automatically generated parameters can be obtained with the [`parameters`](@ref) function
from the [`MonodromyResult`](@ref).

## Options

* `check_startsolutions=true`: If `true`, we do a Newton step for each entry of `sols`for
  checking if it is a valid startsolutions. Solutions which are not valid are sorted out.
* `distance_function=euclidean_distance`: The distance function used for
  [`UniquePoints`](@ref).
* `done_callback=always_false`: A callback to end the computation early. This function takes
  2 arguments. The first one is the new solution `x` and the second one are all current
  solutions (including `x`). Return `true` if the compuation is done.
* `equivalence_classes=true`: This only applies if there is at least one group action
  supplied. We then consider two solutions in the same equivalence class if we can transform
  one to the other by the supplied group actions. We only track one solution per equivalence
  class.
* `identical_tol::Float64=1e-6`: The tolerance with which it is decided whether two
  solutions are identical.
* `group_action=nothing`: A function taking one solution and returning other solutions if
  there is a constructive way to obtain them, e.g. by symmetry.
* `group_actions=nothing`: If there is more than one group action you can use this to chain
  the application of them. For example if you have two group actions `foo` and `bar` you can
  set `group_actions=[foo, bar]`. See [`GroupActions`](@ref) for details regarding the
  application rules.
* `max_loops_no_progress::Int=10`: The maximal number of iterations (i.e. loops generated)
  without any progress.
* `min_solutions`: The minimal number of solutions before a stopping heuristic is applied.
  By default this is half of `target_solutions_count` if applicable otherwise 1.
* `parameter_sampler=independent_normal`: A function taking the parameter `p` and returning
  a new random parameter `q`. By default each entry of the parameter vector is drawn
  independently from the univariate normal distribution.
* `resuse_loops::Symbol=:all`: Strategy to reuse other loops for new found solutions. `:all`
  propagates a new solution through all other loops, `:random` picks a random loop, `:none`
  doesn't reuse a loop.
* `show_progress=true`: Enable a progress meter.
* `strategy`: The strategy used to create loops. If `F` only depends linearly on `p` this
  will be [`Petal`](@ref). Otherwise this will be [`Triangle`](@ref) with weights if `F`
  is a real system.
* `target_solutions_count=nothing`: The computations are stopped if this number of
  solutions is reached.
* `threading = true`: Enable multithreading of the path tracking.
* `timeout=float(typemax(Int))`: The maximal number of *seconds* the computation is allowed
  to run.
"""
function monodromy_solve(
    F::MPPolys;
    parameters = throw(ArgumentError("Necessary to provide `parameters` as a keyword argument.")),
    variable_ordering = nothing,
    seed = randseed(),
    show_progress = true,
    threading::Bool = true,
    kwargs...,
)
    Random.seed!(seed)
    x₀ = p₀ = nothing
    try
        x₀, p₀ = find_start_pair(F, parameters; variable_ordering = variable_ordering)
    catch
        throw(ArgumentError("Cannot compute a start pair (x, p) automatically."))
    end
    monodromy_solve!(
        MonodromySolver(
            F,
            [x₀],
            p₀;
            parameters = parameters,
            variable_ordering = variable_ordering,
            kwargs...,
        ),
        seed;
        show_progress = show_progress,
        threading = threading,
    )
end

function find_start_pair(F::MPPolys, parameters; variable_ordering = nothing)
    if variable_ordering !== nothing
        vars = variable_ordering
    else
        vars = variables(F; parameters = parameters)
    end
    x₀ = randn(ComplexF64, length(vars))
    if F isa Composition
        F = expand(F)
    end
    F_x₀ = [subs(f, vars => x₀) for f in F]
    # let's use the helper function above to get the matrix form
    A, b = linear_system(F_x₀)
    if iszero(b)
        p_end = randn(ComplexF64)
        p₀ = [LA.qr(A[:, 1:end-1], Val(true)) \ ((-p_end) .* A[:, end]); p_end]
    else
        p₀ = LA.qr(A, Val(true)) \ b
    end

    all(
        f -> abs(f(parameters => p₀)) < 1e-8,
        F_x₀,
    ) || throw(ArgumentError("Cannot compute start pair."))

    (x₀, p₀)
end

function monodromy_solve(
    args...;
    seed = randseed(),
    show_progress = true,
    threading::Bool = true,
    kwargs...,
)
    Random.seed!(seed)
    monodromy_solve!(
        MonodromySolver(args...; kwargs...),
        seed;
        show_progress = show_progress,
        threading = threading,
    )
end

function default_strategy(
    F::MPPolyInputs,
    parameters,
    p₀::AbstractVector{TP};
    is_real_system = false,
) where {TC,TP}
    # If F depends only linearly on the parameters a petal is sufficient
    vars = variables(F; parameters = parameters)
    if all(d -> d ≤ 1, maxdegrees(F; parameters = vars))
        Petal()
    # For a real system we should introduce some weights to avoid the discriminant
    elseif is_real_system
        Triangle(useweights = true)
    else
        Triangle(useweights = false)
    end
end
function default_strategy(F::AbstractSystem, parameters, p₀; is_real_system = false)
    # For a real system we should introduce some weights to avoid the discriminant
    if is_real_system
        Triangle(useweights = true)
    else
        Triangle(useweights = false)
    end
end

function numerically_check_real(H::AbstractHomotopy, x)
    y = copy(x)
    Random.randn!(y)
    for i in eachindex(y)
        y[i] = real(y[i]) + 0.0im
    end
    t = rand()
    r = evaluate(H, y, t)
    isapprox(LinearAlgebra.norm(imag.(r), Inf), 0.0, atol = 1e-14)
end

#################
## Actual work ##
#################

"""
    MonodromyJob

A `MonodromyJob` is consisting of a solution id and a loop id.
"""
struct MonodromyJob
    id::Int
    loop_id::Int
end

function monodromy_solve!(
    MS::MonodromySolver,
    seed;
    show_progress::Bool = true,
    threading::Bool = true,
)
    if nsolutions(MS) == 0
        @warn "None of the provided solutions is a valid start solution."
        return MonodromyResult(
            :invalid_startvalue,
            similar(MS.solutions.points, 0),
            MS.parameters,
            MS.statistics,
            MS.options.equivalence_classes,
            seed,
        )
    end
    # solve
    retcode = :not_assigned
    if show_progress
        if !MS.options.equivalence_classes
            desc = "Solutions found:"
        else
            desc = "Classes of solutions (modulo group action) found:"
        end
        progress = ProgressMeter.ProgressUnknown(
            desc;
            delay = 0.5,
            clear_output_ijulia = true,
        )
    else
        progress = nothing
    end

    n_blas_threads = single_thread_blas()

    retcode = monodromy_solve!(MS, seed, progress; threading = threading)

    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)
    finished!(MS.statistics, nsolutions(MS))
    MonodromyResult(
        retcode,
        solutions(MS),
        MS.parameters,
        MS.statistics,
        MS.options.equivalence_classes,
        seed,
    )
end


function monodromy_solve!(
    MS::MonodromySolver,
    seed,
    progress::Union{Nothing,ProgressMeter.ProgressUnknown};
    threading::Bool = true,
)
    Threads.resize_nthreads!(MS.trackers)
    t₀ = time_ns()
    iterations_without_progress = 0 # stopping heuristic
    # intialize job queue
    queue = MonodromyJob.(1:nsolutions(MS), 1)
    thread_queues = [MonodromyJob[] for _ = 1:Threads.nthreads()]

    n = nsolutions(MS)
    retcode = :none
    while n < MS.options.target_solutions_count
        retcode = empty_queue!(
            queue,
            thread_queues,
            MS,
            t₀,
            progress;
            threading = threading,
        )

        if retcode == :done
            retcode = :success
            break
        elseif retcode == :timeout ||
               retcode == :invalid_startvalue || retcode == :interrupt
            break
        end

        for q in thread_queues
            append!(queue, q)
            empty!(q)
        end
        isempty(queue) || continue

        # Iterations heuristic
        n_new = nsolutions(MS)
        if n == n_new
            iterations_without_progress += 1
        else
            iterations_without_progress = 0
            n = n_new
        end
        if iterations_without_progress == MS.options.max_loops_no_progress &&
           n_new ≥ MS.options.min_solutions
            retcode = :heuristic_stop
            break
        end

        regenerate_loop_and_schedule_jobs!(queue, MS)
    end
    update_progress!(progress, nsolutions(MS), MS.statistics; finish = true)

    retcode
end

function empty_queue!(
    queue,
    thread_queues,
    MS::MonodromySolver,
    t₀::UInt,
    progress;
    threading::Bool = true,
)
    if !threading
        ntracked = Ref(0)
        nfailures = Ref(0)
        n_since_update = 0
        retcode = :incomplete
        try
            while !isempty(queue)
                job = pop!(queue)
                status = process!(queue, job, MS, progress, ntracked, nfailures, 1)
                if status == :done || status == :invalid_startvalue || status == :interrupt
                    retcode = status
                    break
                end
                n_since_update += 1

                MS.statistics.ntrackedpaths += ntracked[]
                MS.statistics.ntrackingfailures += nfailures[]
                ntracked[] = nfailures[] = 0
                if n_since_update ≥ 4
                    n_since_update = 0
                    update_progress!(progress, nsolutions(MS), MS.statistics)
                end
                # check timeout
                if (time_ns() - t₀) > MS.options.timeout * 1e9 # convert s to ns
                    retcode = :timeout
                    break
                end
            end
        catch e
            if isa(e, InterruptException)
                return :interrupt
            else
                rethrow()
            end
        end
        n_since_update > 0 && update_progress!(progress, nsolutions(MS), MS.statistics)
        retcode
    else
        stop = Ref(false)
        return_status = Ref(:incomplete)
        start_ntrackedpaths = MS.statistics.ntrackedpaths
        start_ntrackingfailures = MS.statistics.ntrackingfailures
        ntrackedpaths = [Ref(0) for _ = 1:Threads.nthreads()]
        ntrackingfailures = [Ref(0) for _ = 1:Threads.nthreads()]
        Threads.@threads for job in queue
            try
                !stop[] || continue
                tid = Threads.threadid()
                status = process!(
                    thread_queues[tid],
                    job,
                    MS,
                    progress,
                    ntrackedpaths[tid],
                    ntrackingfailures[tid],
                    tid,
                )
                if status == :done
                    return_status[] = :done
                    stop[] = true
                elseif status == :invalid_startvalue
                    return_status[] = :invalid_startvalue
                    stop[] = true
                end
                if tid == 1
                    n_tracked = start_ntrackedpaths + sum(getindex, ntrackedpaths)
                    n_failures = start_ntrackingfailures + sum(getindex, ntrackingfailures)
                    MS.statistics.ntrackedpaths = n_tracked
                    MS.statistics.ntrackingfailures = n_failures
                    update_progress!(progress, nsolutions(MS), MS.statistics)
                end
            # check timeout
                if (time_ns() - t₀) > MS.options.timeout * 1e9 # convert s to ns
                    return_status[] = :invalid_startvalue
                    stop[] = true
                end
                if stop[]
                    break
                end
            catch e
                if isa(e, InterruptException)
                    return_status[] = :interrupt
                    stop[] = true
                    break
                else
                    rethrow()
                end
            end
        end
        empty!(queue)
        n_tracked = start_ntrackedpaths + sum(getindex, ntrackedpaths)
        n_failures = start_ntrackingfailures + sum(getindex, ntrackingfailures)
        MS.statistics.ntrackedpaths = n_tracked
        MS.statistics.ntrackingfailures = n_failures
        update_progress!(progress, nsolutions(MS), MS.statistics)

        return_status[]
    end
end

function process!(
    queue,
    job::MonodromyJob,
    MS::MonodromySolver,
    progress,
    ntrackedpaths::Base.RefValue{Int},
    ntrackingfailures::Base.RefValue{Int},
    tid::Integer,
)
    x = solutions(MS)[job.id]
    loop = MS.loops[job.loop_id]
    retcode = track(MS.trackers[tid], x, loop, ntrackedpaths, ntrackingfailures)
    is_success(retcode) || return :incomplete

    y = pull_back(MS.problem, current_x(MS.trackers[tid]))
    add_and_schedule!(MS, queue, y, job) && return :done

    if MS.options.complex_conjugation
        add_and_schedule!(MS, queue, conj.(y), job) && return :done
    end

    if !MS.options.equivalence_classes
        apply_actions(MS.options.group_actions, y) do yᵢ
            add_and_schedule!(MS, queue, yᵢ, job)
        end
    end
    return :incomplete
end

"""
    add_and_schedule!(MS, queue, y, job)

Add `y` to the current `node` (if it not already exists) and schedule a new job to the
`queue`. Returns `true` if we are done. Otherwise `false`.
"""
function add_and_schedule!(MS::MonodromySolver, queue, y, job::MonodromyJob) where {N,T}
    lock(MS.solutions_lock)
    k = add!(MS.solutions, y, Val(true); tol = MS.options.identical_tol)
    unlock(MS.solutions_lock)
    if k == NOT_FOUND || k == NOT_FOUND_AND_REAL
        # Check if we are done
        isdone(MS.solutions, y, MS.options) && return true
        push!(queue, MonodromyJob(nsolutions(MS), job.loop_id))
        # Schedule also on other loops
        if MS.options.reuse_loops == :random && length(MS.loops) > 1
            r_loop_id = rand(2:length(MS.loops))
            if r_loop_id == job.loop_id
                r_loop_id -= 1
            end
            push!(queue, MonodromyJob(nsolutions(MS), r_loop_id))
        elseif MS.options.reuse_loops == :all
            for r_loop_id = 1:length(MS.loops)
                r_loop_id != job.loop_id || continue
                push!(queue, MonodromyJob(nsolutions(MS), r_loop_id))
            end
        end
    end
    MS.statistics.nreal += (k == NOT_FOUND_AND_REAL)
    false
end

function update_progress!(
    ::Nothing,
    nsolutions,
    statistics::MonodromyStatistics;
    finish = false,
)
    nothing
end
function update_progress!(
    progress,
    nsolutions,
    statistics::MonodromyStatistics;
    finish = false,
)
    ProgressMeter.update!(
        progress,
        nsolutions,
        showvalues = (
            ("# paths tracked", statistics.ntrackedpaths),
            ("# loops generated", statistics.nparametergenerations),
            (
             "# completed loops without change",
             n_completed_loops_without_change(statistics, nsolutions),
            ),
            (
             "# solutions in current loop",
             n_solutions_current_loop(statistics, nsolutions),
            ),
            ("# real solutions", statistics.nreal),
        ),
    )
    if finish
        ProgressMeter.finish!(progress)
    end
    nothing
end

function isdone(solutions::UniquePoints, x, options::MonodromyOptions)
    options.done_callback(x, solutions.points) || length(solutions) ≥ options.target_solutions_count
end

function regenerate_loop_and_schedule_jobs!(queue, MS::MonodromySolver)
    es = MS.loops[1].edges
    loop = MonodromyLoop(
        MS.parameters,
        length(es),
        MS.options,
        weights = !isnothing(es[1].weights),
    )
    push!(MS.loops, loop)
    for id = 1:nsolutions(MS)
        push!(queue, MonodromyJob(id, length(MS.loops)))
    end
    generated_parameters!(MS.statistics, nsolutions(MS))
    nothing
end


##################
## VERIFICATION ##
##################
"""
    verify_solution_completeness(F, res::MonodromyResult;
                                 parameters=..., trace_tol=1e-6, options...)

Verify that the monodromy computation found all solutions by [`monodromy_solve`](@ref).
This uses a multi-projective trace test as described in [^LRS18].
The trace is a numerical value which is 0 if all solutions are found, for this the
`trace_tol` keyword argument is used. The function returns `nothing` if some computation
couldn't be carried out. Otherwise returns a boolean. Note that this function requires the
computation of solutions to another polynomial system using monodromy. This routine can
return `false` although all solutions are found if this additional solution set is not
complete. The `options...` arguments can be everything which is accepted by `solve` and
`monodromy_solve`.

### Example

```
julia> @polyvar x y a b c;

julia> f = x^2+y^2-1;

julia> l = a*x+b*y+c;

julia> res = monodromy_solve([f,l], [-0.6-0.8im, -1.2+0.4im], [1,2,3]; parameters=[a,b,c])
MonodromyResult
==================================
• 2 solutions (0 real)
• return code → heuristic_stop
• 44 tracked paths
• seed → 367230

julia> verify_solution_completeness([f,l], res; parameters=[a,b,c], trace_tol = 1e-8)
[ Info: Compute additional witnesses for completeness...
[ Info: Found 2 additional witnesses
[ Info: Compute trace...
[ Info: Norm of trace: 1.035918995391323e-15
true
```

    verify_solution_completeness(F, S, p; parameters=..., kwargs...)

Verify the solution completeness using the computed solutions `S` to the parameter `p`.

    verify_solution_completeness(TTS, S, W₁₀, p₀::Vector{<:Number}, l₀)

Use the already computeded additional witnesses `W₁₀`. You want to obtain
`TTS`, `W₁₀` and `l₀` as the output from [`solution_completeness_witnesses`](@ref).


[^LRS18]:
    Leykin, Anton, Jose Israel Rodriguez, and Frank Sottile. "Trace test."
    Arnold Mathematical Journal 4.1 (2018): 113-125.
"""
function verify_solution_completeness(F, R::MonodromyResult; kwargs...)
    verify_solution_completeness(F, solutions(R), R.parameters; kwargs...)
end

function verify_solution_completeness(
    F,
    W₀₁::AbstractVector{<:AbstractVector},
    p₀::Vector{<:Number};
    show_progress = true,
    parameters = nothing,
    trace_tol = 1e-6,
    kwargs...,
)
    W₁₀, TTS, l₀ = solution_completeness_witnesses(
        F,
        W₀₁,
        p₀;
        show_progress = show_progress,
        parameters = parameters,
        kwargs...,
    )
    verify_solution_completeness(TTS, W₀₁, W₁₀, p₀, l₀; show_progress = show_progress)
end

function verify_solution_completeness(
    TTS,
    W₀₁::AbstractVector{<:AbstractVector},
    W₁₀::AbstractVector{<:AbstractVector},
    p₀::Vector{<:Number},
    l₀;
    show_progress = true,
    trace_tol = 1e-6,
    kwargs...,
)
    # Combine W₀₁ and W₁₀.
    S = append!([[x; 0.0] for x in W₀₁], W₁₀)
    # To verify that we found all solutions we need move in the pencil
    if show_progress
        @info("Compute trace...")
    end

    trace = track_and_compute_trace(TTS, S, l₀; kwargs...)
    if isnothing(trace)
        return nothing
    else
        show_progress && @info("Norm of trace: $(LinearAlgebra.norm(trace))")
        LinearAlgebra.norm(trace) < trace_tol
    end
end

"""
    solution_completeness_witnesses(F, S, p; parameters=..., kwargs...)

Compute the additional necessary witnesses. Returns a triple `(W₁₀, TTS, l)`
containing the additional witnesses `W₁₀`, a trace test system `TTS` and
the parameters `l` for `TTS`.
"""
function solution_completeness_witnesses(
    F,
    W₀₁,
    p₀::Vector{<:Number};
    parameters = nothing,
    show_progress = true,
    kwargs...,
)
    # generate another start pair
    q₀ = randn(ComplexF64, length(p₀))
    # Construct the trace test system
    TTS = TraceTestSystem(SPSystem(F; parameters = parameters), p₀, q₀ - p₀)

    y₁ = solution(solve(F, W₀₁[1]; p₁ = p₀, p₀ = q₀, parameters = parameters, kwargs...)[1])
    # construct an affine hyperplane l(x) going through y₀
    l₁ = cis.(2π .* rand(length(y₁)))
    push!(l₁, -sum(l₁ .* y₁))
    # This is numerically sometimes not so nice. Let's move to a truly generic one.
    l₀ = randn(ComplexF64, length(l₁))
    y₀ = solution(solve(TTS, [y₁; 1]; p₁ = l₁, p₀ = l₀)[1])

    if show_progress
        @info("Compute additional witnesses for completeness...")
    end

    R₁₀ = monodromy_solve(TTS, y₀, l₀; max_loops_no_progress = 5, kwargs...)
    best_result = R₁₀
    best_params = l₀
    result_agreed = 0
    for i = 1:10
        k₀ = randn(ComplexF64, length(l₀))
        S_k₀ = solutions(solve(
            TTS,
            solutions(R₁₀);
            start_parameters = l₀,
            target_parameters = k₀,
        ))
        new_result = monodromy_solve(TTS, S_k₀, k₀; max_loops_no_progress = 5)
        if nsolutions(new_result) == nsolutions(best_result)
            result_agreed += 1
        elseif nsolutions(new_result) > nsolutions(best_result)
            best_result = new_result
            best_params = k₀
        end
        if result_agreed > 2
            break
        end
    end

    W₁₀ = solutions(best_result)

    if show_progress
        @info("Found $(length(W₁₀)) additional witnesses")
    end
    W₁₀, TTS, best_params
end


function track_and_compute_trace(TTS::TraceTestSystem, S, l₀; kwargs...)
    for i = 1:3
        TTP = TraceTestPencil(TTS, l₀)
        R₁ = solve(TTP, S, start_parameters = [0.0], target_parameters = [0.1], kwargs...)
        R₂ = solve(TTP, S, start_parameters = [0.0], target_parameters = [-.1], kwargs...)
        if nsolutions(R₁) ≠ length(S) || nsolutions(R₂) ≠ length(S)
            if i == 3
                printstyled("Lost solutions $i times. Abort.\n", color = :red, bold = true)
                return nothing
            end
            printstyled(
                "Lost solutions, need to recompute trace...\n",
                color = :yellow,
                bold = true,
            )
        end
        s₁ = sum(solutions(R₁))
        s = sum(S)
        s₂ = sum(solutions(R₂))
        # Due to floating point errors, we compute the RELATIVE trace
        trace = ((s₁ .- s) .- (s .- s₂)) ./ max.(abs.(s₁ .- s), 1.0)
        return trace
    end
    nothing
end
