export monodromy_solve,
    find_start_pair,
    MonodromyOptions,
    MonodromyResult,
    is_success,
    is_heuristic_stop,
    parameters,
    permutations,
    trace,
    verify_solution_completeness

#####################
# Monodromy Options #
#####################


## Options
"""
    MonodromyOptions(; options...)

Options for [`monodromy_solve`](@ref).
"""
Base.@kwdef struct MonodromyOptions{D,GA<:Union{Nothing,GroupActions}}
    check_startsolutions::Bool = true
    group_actions::GA = nothing
    loop_finished_callback = always_false
    parameter_sampler = independent_normal
    equivalence_classes::Bool = !isnothing(group_actions)
    # stopping heuristic
    trace_test_tol::Float64 = 1e-6
    trace_test::Bool = true
    target_solutions_count::Union{Nothing,Int} = nothing
    timeout::Union{Nothing,Float64} = nothing
    min_solutions::Union{Nothing,Int} = nothing
    max_loops_no_progress::Int = 5
    reuse_loops::Symbol = :all
    permutations::Bool = false
    # unique points options
    distance::D = EuclideanNorm()
    triangle_inequality::Union{Nothing,Bool} = nothing
    unique_points_atol::Union{Nothing,Float64} = nothing
    unique_points_rtol::Union{Nothing,Float64} = nothing
    #
    single_loop_per_start_solution = false
end

"""
    independent_normal(p::AbstractVector)

Sample a vector where each entry is drawn independently from the complex Normal distribution
by calling [`randn(ComplexF64)`](@ref).

    independent_normal(L::LinearSubspace)

Creates a random linear subspace by calling [`rand_subspace`](@ref).
"""
independent_normal(p::AbstractVector{T}) where {T} = randn(ComplexF64, length(p))
independent_normal(L::LinearSubspace) = rand_subspace(ambient_dim(L); dim = dim(L))

#############################
# Loops and Data Structures #
#############################

struct LoopTrackingJob
    id::Int
    loop_id::Int
end



#######################
# Loop data structure
#######################


struct MonodromyLoop{P<:Union{LinearSubspace{ComplexF64},Vector{ComplexF64}}}
    # p -> p₁ -> p₂ -> p
    p::P
    p₀₁::P # halfway
    p₁::P
    p₂::P
end

function MonodromyLoop(base::AbstractVector, options::MonodromyOptions)
    p = convert(Vector{ComplexF64}, base)
    p₁ = convert(Vector{ComplexF64}, options.parameter_sampler(p))
    p₂ = convert(Vector{ComplexF64}, options.parameter_sampler(p))

    MonodromyLoop(p, 0.5 .* (p₁ .- p), p₁, p₂)
end

function MonodromyLoop(base::LinearSubspace, options::MonodromyOptions)
    L = convert(LinearSubspace{ComplexF64}, base)
    # the second linear space is just a translation in order to perform a trace test
    # To still find new solutions quickly we translate the linear space
    # by a larger distance
    v = LA.rmul!(LA.normalize!(randn(ComplexF64, codim(L))), 5)
    L₀₁ = translate(L, v, Extrinsic)
    L₁ = translate(L₀₁, v, Extrinsic)
    L₂ = convert(LinearSubspace{ComplexF64}, options.parameter_sampler(L))

    MonodromyLoop(L, L₀₁, L₁, L₂)
end

##########################
## Monodromy Statistics ##
##########################

Base.@kwdef mutable struct MonodromyStatistics
    tracked_loops::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
    tracking_failures::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
    generated_loops::Threads.Atomic{Int} = Threads.Atomic{Int}(0)
    solutions::Vector{Int} = Int[]
    permutations::Vector{Vector{Int}} = Vector{Int}[]
end

MonodromyStatistics(nsolutions::Int) =
    MonodromyStatistics(solutions_development = [nsolutions])


function Base.show(io::IO, S::MonodromyStatistics)
    println(io, "MonodromyStatistics")
    print_fieldname(io, S, :tracked_loops)
    print_fieldname(io, S, :tracking_failures)
    print_fieldname(io, S, :solutions)
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::MonodromyStatistics) = S

# update routines
function loop_tracked!(stats::MonodromyStatistics)
    Threads.atomic_add!(stats.tracked_loops, 1)
    stats
end
function loop_failed!(stats::MonodromyStatistics)
    Threads.atomic_add!(stats.tracking_failures, 1)
    stats
end
function loop_finished!(stats, nsolutions)
    push!(stats.solutions, nsolutions)
end

function add_permutation!(stats::MonodromyStatistics, loop_id, start_id, end_id)
    perms = stats.permutations[loop_id]
    # ensure that perms has correct length
    while length(perms) < start_id
        push!(perms, 0)
    end
    perms[start_id] = end_id
    stats
end

function solutions_current_loop(stats::MonodromyStatistics, nsolutions)
    nsolutions - stats.solutions[end]
end
function solutions_last_loop(stats::MonodromyStatistics)
    length(stats.solutions) > 1 ? stats.solutions[end] - stats.solutions[end-1] : 0
end

function loops_no_change(stats::MonodromyStatistics, nsolutions)
    k = 0
    for i = length(stats.solutions):-1:1
        stats.solutions[i] == nsolutions || break
        k += 1
    end
    return max(k - 1, 0)
end

@noinline function make_showvalues(stats::MonodromyStatistics; queued::Int, solutions::Int)
    [
        ("tracked loops (queued)", "$(stats.tracked_loops[]) ($queued)"),
        (
            "solutions in current (last) loop",
            "$(solutions_current_loop(stats, solutions)) ($(solutions_last_loop(stats)))",
        ),
        (
            "generated loops (no change)",
            "$(stats.generated_loops[]) ($(loops_no_change(stats, solutions)))",
        ),
    ]
end



############
## Result ##
############

"""
    MonodromyResult

The monodromy result contains the result of the [`monodromy_solve`](@ref) computation.
"""
struct MonodromyResult{P,LP}
    returncode::Symbol
    results::Vector{PathResult}
    parameters::P
    loops::Vector{MonodromyLoop{LP}}
    statistics::MonodromyStatistics
    equivalence_classes::Bool
    seed::UInt32
    trace::Union{Nothing,Float64}
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MonodromyResult) = x
function Base.show(io::IO, result::MonodromyResult{T}) where {T}
    println(io, "MonodromyResult")
    println(io, "="^length("MonodromyResult"))
    println(io, "• return_code → :$(result.returncode)")
    if result.equivalence_classes
        println(io, "• $(nsolutions(result)) classes of solutions (modulo group action)")
    else
        println(io, "• $(nsolutions(result)) solutions")
    end
    println(io, "• $(result.statistics.tracked_loops[]) tracked loops")
    print(io, "• random_seed → ", sprint(show, result.seed))
    if !isnothing(result.trace)
        print(io, "\n• trace → ", sprint(show, result.trace))
    end
end

TreeViews.hastreeview(::MonodromyResult) = true
TreeViews.numberofnodes(::MonodromyResult) = 7
TreeViews.treelabel(io::IO, x::MonodromyResult, ::MIME"application/prs.juno.inline") =
    print(
        io,
        "<span class=\"syntax--support syntax--type syntax--julia\">MonodromyResult</span>",
    )
#
function TreeViews.nodelabel(
    io::IO,
    x::MonodromyResult,
    i::Int,
    ::MIME"application/prs.juno.inline",
)
    if i == 1
        if x.equivalence_classes
            print(io, "$(nsolutions(x)) solutions (modulo group action)")
        else
            print(io, "$(nsolutions(x)) solutions")
        end
    elseif i == 2
        print(io, "Return code")
    elseif i == 3
        print(io, "Tracked loops")
    elseif i == 4
        print(io, "Statistics")
    elseif i == 5
        print(io, "Parameters")
    elseif i == 6
        print(io, "Seed")
    elseif i == 7
        print(io, "Trace")
    end
end
function TreeViews.treenode(r::MonodromyResult, i::Integer)
    if i == 1
        return r.results
    elseif i == 2
        return r.returncode
    elseif i == 3
        return r.statistics.tracked_loops[]
    elseif i == 4
        return r.statistics
    elseif i == 5
        return r.parameters
    elseif i == 6
        return r.seed
    elseif i == 7
        return something(r.trace, missing)
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
    solutions(result::MonodromyResult)

Return all solutions.
"""
solutions(r::MonodromyResult) = solutions(r.results)

"""
    nsolutions(result::MonodromyResult)

Returns the number solutions of the `result`.
"""
nsolutions(r::MonodromyResult) = length(r.results)

"""
    results(result::MonodromyResult)

Returns the computed [`PathResult`](@ref)s.
"""
results(r::MonodromyResult) = r.results

"""
    nresults(result::MonodromyResult)

Returns the number of results computed.
"""
nresults(r::MonodromyResult) = length(r.results)

"""
    parameters(result::MonodromyResult)

Return the parameters corresponding to the given result `r`.
"""
ModelKit.parameters(r::MonodromyResult) = r.parameters

"""
    seed(result::MonodromyResult)

Return the random seed used for the computations.
"""
seed(r::MonodromyResult) = r.seed

"""
    trace(result::MonodromyResult)

Return the result of the trace test computed during the monodromy.
"""
trace(r::MonodromyResult) = r.trace


"""
    permutations(r::MonodromyResult; reduced=true)

Return the permutations of the solutions that are induced by tracking over the loops.
If `reduced=false`, then all permutations are returned.
If `reduced=true` then permutations without repetitions are returned.

If a solution was not tracked in the loop, then the corresponding entry is 0.

Example: monodromy loop for a varying line that intersects two circles.
```julia
using LinearAlgebra
@var x[1:2] a b c
c1 = (x - [2, 0]) ⋅ (x - [2, 0]) - 1
c2 = (x - [-2, 0]) ⋅ (x - [-2, 0]) - 1
F = [c1 * c2; a * x[1] + b * x[2] - c]
S = monodromy_solve(F, [[1, 0]], [1, 1, 1], parameters = [a, b, c], permutations = true)

permutations(S)
```

will return

```julia
2×2 Array{Int64,2}:
 1  2
 2  1
```

and `permutations(S, reduced = false)` returns

```julia
2×12 Array{Int64,2}:
 1  2  2  1  1  …  1  2  1  1  1
 2  1  1  2  2     2  1  2  2  2
```

"""
function permutations(r::MonodromyResult; reduced::Bool = true)
    π = r.statistics.permutations
    N = nresults(r)

    π = filter(πⱼ -> length(πⱼ) == N, π)
    if reduced
        π = unique(π)
    end

    A = zeros(Int, N, length(π))
    for (j, πⱼ) in enumerate(π), i = 1:N
        A[i, j] = πⱼ[i]
    end

    A
end

#####################
## monodromy solve ##
#####################

mutable struct MonodromySolver{
    Tracker<:EndgameTracker,
    P,
    UP<:UniquePoints,
    MO<:MonodromyOptions,
}
    trackers::Vector{Tracker}
    loops::Vector{MonodromyLoop{P}}
    unique_points::UP
    unique_points_lock::ReentrantLock
    options::MO
    statistics::MonodromyStatistics
    # We save the sums of the solutions for three
    # values to check if the trace is colinear
    # We augment the sums by 1 in each coordinate to
    # make this work for the two variables case
    # (n + 1) × 3 matrix
    trace::Matrix{ComplexF64}
    trace_lock::ReentrantLock
end

function MonodromySolver(
    F::Union{System,AbstractSystem},
    parameters;
    intrinsic = nothing,
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    tracker_options = TrackerOptions(),
    options = MonodromyOptions(),
)
    if !isa(options, MonodromyOptions)
        options = MonodromyOptions(; options...)
    end

    if parameters isa LinearSubspace
        H = linear_subspace_homotopy(
            F,
            parameters,
            parameters;
            compile = compile,
            intrinsic = intrinsic,
        )
        P = LinearSubspace{ComplexF64}
    else
        H = parameter_homotopy(F; generic_parameters = parameters, compile = compile)
        P = Vector{ComplexF64}
    end
    egtracker = EndgameTracker(
        H;
        tracker_options = tracker_options,
        options = EndgameOptions(endgame_start = 0.0),
    )
    trackers = [egtracker]
    group_actions = options.equivalence_classes ? options.group_actions : nothing
    x₀ = zeros(ComplexF64, size(F, 2))
    #TODO: This is wrong for IntrinsicSubspaceHomotopy of projective varieties
    if !isnothing(group_actions) && (egtracker.tracker.homotopy isa AffineChartHomotopy)
        group_actions = let H = egtracker.tracker.homotopy, group_actions = group_actions
            (cb, s) -> apply_actions(v -> cb(set_solution!(v, H, v, 0)), group_actions, s)
        end
    end

    unique_points = UniquePoints(
        x₀,
        1;
        distance = options.distance,
        group_actions = group_actions,
        triangle_inequality = options.triangle_inequality,
    )

    trace = zeros(ComplexF64, length(x₀) + 1, 3)
    MonodromySolver(
        trackers,
        MonodromyLoop{P}[],
        unique_points,
        ReentrantLock(),
        options,
        MonodromyStatistics(),
        trace,
        ReentrantLock(),
    )
end

function add_loop!(MS::MonodromySolver, p)
    push!(MS.loops, MonodromyLoop(p, MS.options))
    Threads.atomic_add!(MS.statistics.generated_loops, 1)
    if MS.options.permutations
        push!(MS.statistics.permutations, zeros(Int, length(MS.unique_points)))
    end
    MS
end
loop(MS::MonodromySolver, i) = MS.loops[i]
nloops(MS::MonodromySolver) = length(MS.loops)

function reset_loops!(MS::MonodromySolver)
    empty!(MS.loops)
    MS
end

function reset_trace!(MS::MonodromySolver)
    MS.trace .= 0
    MS.trace[end, :] .= 1
    MS
end

# Check if trace is colinear
function trace_colinearity(MS)
    σ = LA.svdvals(MS.trace)
    σ[3] / σ[1]
end

"""
    monodromy_solve(F, [sols, p]; options..., tracker_options = TrackerOptions())

Solve a polynomial system `F(x;p)` with specified parameters and initial solutions `sols`
by monodromy techniques. This makes loops in the parameter space of `F` to find new solutions.
If `F` the parameters `p` only occur *linearly* in `F` it is eventually possible to compute
a *start pair* ``(x₀, p₀)`` automatically. In this case `sols` and `p` can be omitted and
the automatically generated parameters can be obtained with the [`parameters`](@ref) function
from the [`MonodromyResult`](@ref).

    monodromy_solve(F, [sols, L]; dim, codim, intrinsic = nothing, options...,
                                  tracker_options = TrackerOptions())

Solve the polynomial system `[F(x); L(x)] = 0` where `L` is a `[LinearSubspace`](@ref).
If `sols` and `L` are not provided it is necesary to provide `dim` or `codim` where `(co)dim`
is the expected (co)dimension of a component of `V(F)`.
See also [`linear_subspace_homotopy`](@ref) for the `intrinsic` option.

## Options

* `catch_interrupt = true`: If true catches interruptions (e.g. issued by pressing Ctrl-C)
  and returns the partial result.
* `check_startsolutions = true`: If `true`, we do a Newton step for each entry of `sols`for
  checking if it is a valid startsolutions. Solutions which are not valid are sorted out.
* `compile = $(COMPILE_DEFAULT[])`: If `true` then a `System` (resp. `Homotopy`) is compiled
  to a straight line program ([`CompiledSystem`](@ref) resp. [`CompiledHomotopy`](@ref))
  for evaluation. This induces a compilation overhead. If `false` then the generated program
  is only interpreted ([`InterpretedSystem`](@ref) resp. [`InterpretedHomotopy`](@ref)).
  This is slower than the compiled version, but does not introduce compilation overhead.
* `distance = EuclideanNorm()`: The distance function used for [`UniquePoints`](@ref).
* `loop_finished_callback = always_false`: A callback to end the computation. This function is
  called with all current [`PathResult`](@ref)s after a loop is exhausted.
  2 arguments. Return `true` if the compuation should be stopped.
* `equivalence_classes=true`: This only applies if there is at least one group action
  supplied. We then consider two solutions in the same equivalence class if we can transform
  one to the other by the supplied group actions. We only track one solution per equivalence
  class.
* `group_action = nothing`: A function taking one solution and returning other solutions if
  there is a constructive way to obtain them, e.g. by symmetry.
* `group_actions = nothing`: If there is more than one group action you can use this to chain
  the application of them. For example if you have two group actions `foo` and `bar` you can
  set `group_actions=[foo, bar]`. See [`GroupActions`](@ref) for details regarding the
  application rules.
* `max_loops_no_progress = 5`: The maximal number of iterations (i.e. loops generated)
  without any progress.
* `min_solutions`: The minimal number of solutions before a stopping heuristic is applied.
  By default no lower limit is enforced.
* `parameter_sampler = independent_normal`: A function taking the parameter `p` and
  returning a new random parameter `q`. By default each entry of the parameter vector
  is drawn independently from Normal distribution.
* `permutations = false`: Whether to keep track of the permutations induced by the loops.
* `resuse_loops = :all`: Strategy to reuse other loops for new found solutions. `:all`
  propagates a new solution through all other loops, `:random` picks a random loop, `:none`
  doesn't reuse a loop.
* `target_solutions_count`: The computation is stopped if this number of solutions is
  reached.
* `threading = true`: Enable multithreading of the path tracking.
  Careful: Some CPUs hang when using multiple threads. To avoid this run Julia with 1 
  interactive thread for the REPL and `n` threads for other tasks (e.g., `julia -t 8,1` for `n=8`).
* `timeout`: The maximal number of *seconds* the computation is allowed to run.
* `trace_test = true`: If `true` a trace test is performed to check whether all solutions
  are found. This is only applicable if monodromy is performed with a linear subspace.
  See also [`trace_test`](@ref).
* `trace_test_tol = 1e-10`: The tolerance for the trace test to be successfull.
  The trace is divided by the number of solutions before compared to the trace_test_tol.
* `unique_points_rtol`: the relative tolerance for [`unique_points`](@ref).
* `unique_points_atol`: the absolute tolerance for [`unique_points`](@ref).
"""
function monodromy_solve(F::Vector{<:ModelKit.Expression}, args...; parameters, kwargs...)
    monodromy_solve(System(F; parameters = parameters), args...; kwargs...)
end
function monodromy_solve(
    F::AbstractVector{<:MP.AbstractPolynomial},
    args...;
    parameters,
    variables = setdiff(MP.variables(F), parameters),
    variable_ordering = variables,
    kwargs...,
)
    sys = System(F, variables = variable_ordering, parameters = parameters)
    monodromy_solve(sys, args...; kwargs...)
end
function monodromy_solve(
    F::Union{System,AbstractSystem},
    args...;
    seed = rand(UInt32),
    tracker_options = TrackerOptions(),
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    catch_interrupt::Bool = true,
    dim = nothing,
    codim = nothing,
    intrinsic = nothing,
    # monodromy options
    options = nothing,
    group_action = nothing,
    group_actions = isnothing(group_action) ? nothing : GroupActions(group_action),
    kwargs...,
)
    if isnothing(options)
        if !(group_actions isa GroupActions) && !isnothing(group_actions)
            group_actions = GroupActions(group_actions)
        end
        options = MonodromyOptions(; group_actions = group_actions, kwargs...)
    end

    if !isnothing(seed)
        Random.seed!(seed)
    end
    if length(args) == 0
        start_pair = find_start_pair(fixed(F; compile = compile))
        if isnothing(start_pair)
            error(
                "Cannot compute a start pair (x, p) using `find_start_pair(F)`." *
                " You need to explicitly pass a start pair.",
            )
        end
        x, p = start_pair
        S = [x]
        # if we have no parameters, then we intersect with a linear subspace.
        # For this we need to know the intended dimension or codimension.
        # We don't guess dim to catch the case that the user just forgot to pass
        # parameters
        if isnothing(p)
            if isnothing(dim) && isnothing(codim)
                error(
                    "Given system doesn't have any parameters. If you intended to intersect " *
                    "with a linear subspace it is necessary to provide a " *
                    "dimension (`dim`) or codimension (`codim`) of the component of interest.",
                )
            else
                projective = is_homogeneous(System(F))
                if !isnothing(codim)
                    codim += projective
                end
                p = rand_subspace(
                    x;
                    dim = codim,
                    codim = dim,
                    affine = !is_homogeneous(System(F)),
                )
            end
        end
    else
        S, p = args
    end
    MS = MonodromySolver(
        F,
        p;
        compile = compile,
        options = options,
        tracker_options = tracker_options,
        intrinsic = intrinsic,
    )

    monodromy_solve(
        MS,
        S,
        p,
        seed;
        show_progress = show_progress,
        threading = threading,
        catch_interrupt = catch_interrupt,
    )
end

"""
    find_start_pair(F; max_tries = 100, atol = 0.0, rtol = 1e-12)

Try to find a pair `(x,p)` for the system `F` such that `F(x,p) = 0` by randomly sampling
a pair `(x₀, p₀)` and performing Newton's method in variable *and* parameter space.
"""
function find_start_pair(
    F::System;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    max_tries = 1_000,
    kwargs...,
)
    if nparameters(F) != 0
        start_pair = find_start_pair_linear_in_params(F; kwargs...)
        if !isnothing(start_pair)
            return start_pair
        end
    end
    find_start_pair(fixed(F; compile = compile); max_tries = max_tries, kwargs...)
end
function find_start_pair(
    F::AbstractSystem;
    max_tries::Int = 1_000,
    atol::Float64 = 0.0,
    rtol::Float64 = 1e-12,
)
    if nparameters(F) == 0
        SF = F
    else
        SF = StartPairSystem(F)
    end
    m = nparameters(F)
    n = nvariables(F)
    c = NewtonCache(SF)
    for i = 1:max_tries
        xp₀ = randn(ComplexF64, size(SF, 2))
        res = newton(SF, xp₀, nothing, InfNorm(), c, atol = atol, rtol = rtol)
        if is_success(res)
            x, p = res.x[1:n], res.x[n+1:end]
            res2 = newton(F, x, p, InfNorm())
            if is_success(res2)
                return res2.x, isempty(p) ? nothing : p
            end
        end
    end
    nothing
end
function find_start_pair_linear_in_params(F; kwargs...)
    x₀ = randn(ComplexF64, nvariables(F))
    F_x₀ = F(x₀, parameters(F))
    all(!ModelKit.is_number, F_x₀) || return nothing

    Ab = construct_linear_system(F_x₀, parameters(F))
    if !isnothing(Ab)
        A, b = Ab
        m, n = size(A)
        m ≤ n || return nothing
        # if b is 0 then we compute the a basis vector of the nullspace
        if iszero(b)
            # don't want to have 0 as a parameter
            m == n && return nothing
            N = LA.nullspace(A)
            p₀ = N * randn(ComplexF64, size(N, 2))
        else
            p₀ = LA.qr(A, Val(true)) \ b
        end
        return x₀, p₀
    else
        return nothing
    end
end

function construct_linear_system(F, x)
    A = zeros(ComplexF64, length(F), length(x))
    b = zeros(ComplexF64, length(F))

    try
        for (i, f) in enumerate(F)
            M, coeffs = exponents_coefficients(f, x)
            for (j, exp) in enumerate(eachcol(M))
                var_found = false
                for (k, d) in enumerate(exp)
                    if d < 0 || d > 1 || (d == 1 && var_found)
                        return nothing
                    end
                    if d == 1
                        var_found = true
                        A[i, k] = coeffs[j]
                    end
                end
                if !var_found
                    b[i] = -coeffs[j]
                end
            end
        end
    catch err
        if err isa ModelKit.PolynomialError
            return nothing
        else
            rethrow(err)
        end
    end
    A, b
end

monodromy_solve(MS::MonodromySolver, x::AbstractVector{<:Number}, p, seed; kwargs...) =
    monodromy_solve(MS, [x], p, seed; kwargs...)

function monodromy_solve(
    MS::MonodromySolver,
    X::AbstractVector{<:AbstractVector{<:Number}},
    p,
    seed;
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
)
    if !show_progress
        progress = nothing
    else
        if MS.options.equivalence_classes
            desc = "Solutions (modulo group action) found:"
        else
            desc = "Solutions found:"
        end
        progress = ProgressMeter.ProgressUnknown(; dt = 0.4, desc = desc, output = stdout)
        progress.tlast += 0.3
    end
    MS.statistics = MonodromyStatistics()
    empty!(MS.unique_points)
    reset_trace!(MS)
    reset_loops!(MS)
    results = check_start_solutions(MS, X, p)
    if isempty(results)
        @warn "None of the provided solutions is a valid start solution (Newton's method did not converge)."
        retcode = :invalid_startvalue
    else
        retcode = :default
        try
            # convert to complex base
            if p isa LinearSubspace
                cp = convert(LinearSubspace{ComplexF64}, p)
            else
                cp = convert(Vector{ComplexF64}, p)
            end
            if threading
                retcode = threaded_monodromy_solve!(results, MS, cp, seed, progress)
            else
                retcode = serial_monodromy_solve!(results, MS, cp, seed, progress)
            end
        catch e
            if !catch_interrupt || !(
                isa(e, InterruptException) ||
                (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
            )
                rethrow(e)
            end
            retcode = :interrupted
        end
    end

    MonodromyResult(
        retcode,
        results,
        p,
        MS.loops,
        MS.statistics,
        MS.options.equivalence_classes,
        seed,
        p isa LinearSubspace ? trace_colinearity(MS) : nothing,
    )
end

function check_start_solutions(MS::MonodromySolver, X, p)
    tracker = MS.trackers[1]
    parameters!(tracker, p, p)
    results = PathResult[]
    for x in X
        res = track(tracker, x)
        if is_success(res)
            _, added = add!(MS, res, length(results) + 1)
            if added
                push!(results, res)
            end
        end
    end

    results
end

function add!(MS::MonodromySolver, res::PathResult, id)
    if isnothing(MS.options.unique_points_rtol)
        rtol = uniqueness_rtol(res)
    else
        rtol = MS.options.unique_points_rtol
    end
    if isnothing(MS.options.unique_points_atol)
        atol = 1e-14
    else
        atol = MS.options.unique_points_atol
    end

    add!(MS.unique_points, solution(res), id; atol = atol, rtol = rtol)
end

function serial_monodromy_solve!(
    results::Vector{PathResult},
    MS::MonodromySolver,
    p,
    seed,
    progress,
)
    queue = LoopTrackingJob[]
    tracker = MS.trackers[1]
    t₀ = time()
    retcode = :in_progress
    stats = MS.statistics
    while true
        loop_finished!(stats, length(results))

        if MS.options.loop_finished_callback(results)
            retcode = :terminated_callback
            @goto _return
        end
        if p isa LinearSubspace &&
           nloops(MS) > 0 &&
           MS.options.trace_test &&
           trace_colinearity(MS) < MS.options.trace_test_tol
            retcode = :success
            @goto _return
        end
        if isnothing(MS.options.target_solutions_count) &&
           loops_no_change(stats, length(results)) ≥ MS.options.max_loops_no_progress
            retcode = :heuristic_stop
            @goto _return
        end
        if length(results) == something(MS.options.target_solutions_count, -1)
            retcode = :success
            @goto _return
        end

        if nloops(MS) > 0 && MS.options.single_loop_per_start_solution
            retcode = :success
            @goto _return
        end

        add_loop!(MS, p)
        reset_trace!(MS)
        # schedule all jobs
        new_loop_id = nloops(MS)
        for i = 1:length(results)
            push!(queue, LoopTrackingJob(i, new_loop_id))
        end

        while !isempty(queue)
            job = popfirst!(queue)
            res = track(
                tracker,
                results[job.id],
                loop(MS, job.loop_id),
                MS.trace,
                MS.trace_lock;
                collect_trace = MS.options.trace_test && nloops(MS) == job.loop_id,
            )
            if !isnothing(res)
                loop_tracked!(stats)

                # 1) check whether solutions already exists
                id, got_added = add!(MS, res, length(results) + 1)

                if MS.options.permutations
                    add_permutation!(stats, job.loop_id, job.id, id)
                end

                got_added || @goto _update
                # 2) doesn't exist, so add to results
                push!(results, res)

                # 3) schedule on same loop again
                if !MS.options.single_loop_per_start_solution
                    push!(queue, LoopTrackingJob(id, job.loop_id))
                end

                # 4) schedule on other loops
                if MS.options.reuse_loops == :all
                    for k = 1:nloops(MS)
                        k != job.loop_id || continue
                        push!(queue, LoopTrackingJob(id, k))
                    end
                elseif MS.options.reuse_loops == :random && nloops(MS) ≥ 2
                    k = rand(2:nloops(MS))
                    if k ≤ job.loop_id
                        k -= 1
                    end
                    push!(queue, LoopTrackingJob(id, k))
                end
            else
                loop_failed!(stats)
                if MS.options.permutations
                    add_permutation!(stats, job.loop_id, job.id, 0)
                end
            end
            # Update progress
            @label _update
            update_progress!(
                progress,
                stats;
                solutions = length(results),
                queued = length(queue),
            )

            if length(results) == something(MS.options.target_solutions_count, -1) &&
               # only terminate after a completed loop to ensure that we collect proper
               # permutation informations
               !MS.options.permutations
                retcode = :success
                @goto _return
            elseif !isnothing(MS.options.timeout) &&
                   time() - t₀ > (MS.options.timeout::Float64)
                retcode = :timeout
                @goto _return
            end
        end
    end

    @label _return

    update_progress!(
        progress,
        stats;
        finish = true,
        solutions = length(results),
        queued = length(queue),
    )

    return retcode
end

function threaded_monodromy_solve!(
    results::Vector{PathResult},
    MS::MonodromySolver,
    p,
    seed,
    progress,
)
    # Shared queue of work
    queue = Channel{LoopTrackingJob}(Inf)

    # Ensure we have one per-thread tracker
    tracker = MS.trackers[1]
    ntrackers = length(MS.trackers)
    nthr = Threads.nthreads()
    resize!(MS.trackers, nthr)
    for i = ntrackers+1:nthr
        MS.trackers[i] = deepcopy(tracker)
    end

    data_lock = ReentrantLock()
    t₀ = time()
    retcode = :in_progress
    stats = MS.statistics
    notify_lock = ReentrantLock()
    interrupted = Ref(false)
    queued = Ref(0)

    progress_lock = ReentrantLock()
    terminated = Threads.Atomic{Bool}(false)

    # Single barrier: spawn workers and controller inside
    try
        Threads.@sync begin
            # Spawn one worker per tracker, but capture a fresh copy per iteration
            for tracker in MS.trackers
                let local_tracker = tracker
                    Threads.@spawn begin
                        for job in queue
                            if terminated[]
                                break
                            end

                            r = track(
                                local_tracker,
                                results[job.id],
                                loop(MS, job.loop_id),
                                MS.trace,
                                MS.trace_lock;
                                collect_trace = MS.options.trace_test && nloops(MS) == job.loop_id
                            )

                            if !isnothing(r)
                                loop_tracked!(stats)
                                lock(data_lock) do
                                    id, added = add!(MS, r, length(results) + 1)
                                    if MS.options.permutations
                                        add_permutation!(stats, job.loop_id, job.id, id)
                                    end
                                    if added
                                        push!(results, r)
                                    end
                                end
                            # Schedule new work as needed
                            else
                                loop_failed!(stats)
                                if MS.options.permutations
                                    add_permutation!(stats, job.loop_id, job.id, 0)
                                end
                            end

                            @lock progress_lock begin
                                update_progress!(progress, stats; solutions = length(results), queued = Base.n_avail(queue))
                            end

                            if length(results) == something(MS.options.target_solutions_count, -1) &&
                            !MS.options.permutations
                                retcode = :success
                                Base.@lock notify_lock begin
                                    interrupted[] = true
                                end
                                try
                                    close(queue)
                                catch e
                                # ignore if already closed
                                end
                            elseif !isnothing(MS.options.timeout) &&
                                time() - t₀ > (MS.options.timeout::Float64)
                                retcode = :timeout
                                Base.@lock notify_lock begin
                                    interrupted[] = true
                                end
                                try
                                    close(queue)
                                catch e
                                # ignore if already closed
                                end
                        end
                        end
                    end
                end
            end

            # Controller thread block - adds loops and schedules work
            Threads.@spawn begin
                Base.@lock notify_lock begin
                    while true
                        if interrupted[]
                            break
                        end
                        loop_finished!(stats, length(results))

                        if MS.options.loop_finished_callback(results)
                            retcode = :terminated_callback
                            break
                        end
                        if loops_no_change(stats, length(results)) ≥ MS.options.max_loops_no_progress
                            retcode = :heuristic_stop
                            break
                        end
                        if length(results) == something(MS.options.target_solutions_count, -1)
                            retcode = :success
                            break
                        end
                        if p isa LinearSubspace &&
                        nloops(MS) > 0 &&
                        MS.options.trace_test &&
                        trace_colinearity(MS) < MS.options.trace_test_tol
                            retcode = :success
                            break
                        end

                        add_loop!(MS, p)
                        reset_trace!(MS)
                        new_loop_id = nloops(MS)
                        for i = 1:length(results)
                            push!(queue, LoopTrackingJob(i, new_loop_id))
                        end

        
                        if (MS.options.single_loop_per_start_solution)
                            retcode = :success
                            break
                        end
                        if retcode != :in_progress
                            break
                        end
                        sleep(0.001)
                    end
                end
            end
        end
    catch e
        close(queue)
        interrupted[] = true
        rethrow(e)
    finally
        queued[] = Base.n_avail(queue)
        close(queue)
    end

    update_progress!(
        progress,
        stats;
        finish = true,
        solutions = length(results),
        queued = queued[],
    )

    return retcode
end

function uniqueness_rtol(res::PathResult)
    # we should have a region of uniqueness of 0.5inv(ω) * norm(x)
    # the norm(x) factors comes from the fact that we computed with a
    # weighted norm.
    # To be more pessimistic we only consider 0.25inv(ω)^2
    # and require at most sqrt(res.accuracy) and at least 1e-14.
    clamp(0.25 * inv(res.ω)^2, 1e-14, sqrt(res.accuracy))
end

"""
    track(tracker, x, edge::LoopEdge, loop::MonodromyLoop, stats::MonodromyStatistics)

Track `x` along the edge `edge` in the loop `loop` using `tracker`. Record statistics
in `stats`.
"""
function track(
    egtracker::EndgameTracker,
    r::PathResult,
    loop::MonodromyLoop,
    trace,
    trace_lock;
    collect_trace::Bool = false,
)
    tracker = egtracker.tracker

    if loop.p isa LinearSubspace && collect_trace
        parameters!(tracker, loop.p, loop.p₀₁)
        x₀ = solution(r)
        retcode =
            track!(tracker, x₀; ω = r.ω, μ = r.μ, extended_precision = r.extended_precision)
        x₀₁ = solution(tracker)
        if !is_success(retcode)
            return nothing
        end

        parameters!(tracker, loop.p₀₁, loop.p₁)
        retcode = track!(
            tracker,
            x₀₁;
            ω = tracker.state.ω,
            μ = tracker.state.μ,
            extended_precision = tracker.state.extended_prec,
        )
        x₁ = solution(tracker)
        if !is_success(retcode)
            return nothing
        end

        Base.@lock trace_lock begin
            for i = 1:length(x₀)
                trace[i, 1] += x₀[i]
                trace[i, 2] += x₀₁[i]
                trace[i, 3] += x₁[i]
            end
        end
    else
        parameters!(tracker, loop.p, loop.p₁)
        retcode = track!(
            tracker,
            solution(r);
            ω = r.ω,
            μ = r.μ,
            extended_precision = r.extended_precision,
        )
        if !is_success(retcode)
            return nothing
        end
        x₁ = solution(tracker)
    end

    parameters!(tracker, loop.p₁, loop.p₂)
    retcode = track!(
        tracker,
        x₁;
        ω = tracker.state.ω,
        μ = tracker.state.μ,
        extended_precision = tracker.state.extended_prec,
    )
    if !is_success(retcode)
        return nothing
    end

    x₂ = solution(tracker)
    parameters!(tracker, loop.p₂, loop.p)
    result = track(
        egtracker,
        x₂;
        ω = tracker.state.ω,
        μ = tracker.state.μ,
        extended_precision = tracker.state.extended_prec,
    )

    if is_success(result)
        result
    else
        nothing
    end
end

update_progress!(progress::Nothing, stats::MonodromyStatistics; kwargs...) = nothing
function update_progress!(
    progress,
    stats::MonodromyStatistics;
    queued::Int,
    solutions::Int,
    finish::Bool = false,
)
    if finish
        ProgressMeter.update!(progress, solutions)
        showvalues = make_showvalues(stats; queued = queued, solutions = solutions)
        ProgressMeter.finish!(progress; showvalues = showvalues)
    elseif time() > progress.tlast + progress.dt
        showvalues = make_showvalues(stats; queued = queued, solutions = solutions)
        ProgressMeter.update!(progress, solutions, showvalues = showvalues)
    end
    # yield()
    nothing
end

## Parameter homotopy constructor
solve(F, R::MonodromyResult; target_parameters, kwargs...) = solve(
    F,
    solutions(R);
    start_parameters = parameters(R),
    target_parameters = target_parameters,
    kwargs...,
)

##################
## VERIFICATION ##
##################
"""
    verify_solution_completeness(F::System, monodromy_result; options...)
    verify_solution_completeness(F::System, solutions, parameters;
        trace_tol = 1e-14,
        show_progress = true,
        compile = COMPILE_DEFAULT[],
        monodromy_options = (compile = compile,),
        parameter_homotopy_options = (compile = compile,),
    )

Verify that a monodromy computation found all solutions by [`monodromy_solve`](@ref).
This uses the trace test described in [^dCR17] and [^LRS18].
The trace is a numerical value which is 0 if all solutions are found, for this the
`trace_tol` keyword argument is used. The function returns `nothing` if some computation
couldn't be carried out. Otherwise returns a boolean. Note that this function requires the
computation of solutions to another polynomial system using monodromy. This routine can
return `false` although all solutions are found if this additional solution set is not
complete.

### Example
```julia
@var x y a b c;
f = x^2+y^2-1;
l = a*x+b*y+c;
sys = System([f, l]; parameters = [a, b, c]);
mres = monodromy_solve(sys, [-0.6-0.8im, -1.2+0.4im], [1,2,3]);
show(mres);
verify_solution_completeness(sys, mres)
```
```
MonodromyResult
==================================
• 2 solutions (0 real)
• return code → heuristic_stop
• 44 tracked paths
• seed → 367230
julia> verify_solution_completeness(sys, mres)
[ Info: Certify provided solutions...
[ Info: Got 2 dinstinct solutions.
[ Info: Compute additional witnesses for completeness...
┌ Info: MonodromyResult
│ ===============
│ • return_code → :heuristic_stop
│ • 4 solutions
│ • 28 tracked loops
└ • random_seed → 0x21e7406a
[ Info: Certify additional witnesses...
[ Info: Computed 2 additional witnesses
[ Info: Compute trace using two parameter homotopies...
[ Info: Norm of trace: 9.33238819760471e-17
true
```

[^dCR17]:
    del Campo, Abraham Martín, and Jose Israel Rodriguez.
    "Critical points via monodromy and local methods."
    Journal of Symbolic Computation 79 (2017): 559-574.

[^LRS18]:
    Leykin, Anton, Jose Israel Rodriguez, and Frank Sottile. "Trace test."
    Arnold Mathematical Journal 4.1 (2018): 113-125.
"""
verify_solution_completeness(F::System, mres::MonodromyResult; kwargs...) =
    verify_solution_completeness(F, solutions(mres), parameters(mres); kwargs...)
function verify_solution_completeness(
    F::System,
    sols::AbstractVector{<:AbstractVector},
    q::AbstractVector;
    show_progress = true,
    compile = COMPILE_DEFAULT[],
    monodromy_options = (compile = compile,),
    parameter_homotopy_options = (compile = compile,),
    trace_tol = 1e-14,
)
    n = nvariables(F)
    m = nparameters(F)
    @unique_var t v[1:m] a[1:n] λ

    x = variables(F)
    p = parameters(F)

    verify_system = System(
        [F(x, p + λ * v); (sum(a .* x) - 1) * λ + t],
        variables = [x; λ],
        parameters = [t; p; v; a],
    )
    # perform monodromy computation to compute the additional witnesses
    if show_progress
        @info "Compute additional witnesses for completeness check..."
    end
    # use the verify_system but enforce t = 0 and start with λ ≠ 0
    # this way we stay on a different irreducible component

    # Let's compute some start solutions for some parameters
    Y, base_params = let
        # We sample some random parameters to set v = qq - q
        qq = randn(ComplexF64, length(q))
        # Its good to have more than 1 start solution, we can construct n
        # by performing a parameter homotopy to qq
        # and then  computing an `a` such that n solutions are on the linear space  a⋅x-1
        qq_res = solve(
            F,
            sols[1:min(n, length(sols))];
            start_parameters = q,
            target_parameters = qq,
            show_progress = show_progress,
            parameter_homotopy_options...,
        )
        a0 = reduce(vcat, transpose.(solutions(qq_res))) \ ones(nsolutions(qq_res))
        map(s -> [s; 1], solutions(qq_res)), [q; qq - q; a0]
    end

    # now use monodromy to find more solutions
    additional_mres = monodromy_solve(
        verify_system,
        Y,
        [0.0; base_params];
        parameter_sampler = p -> [0; randn(ComplexF64, length(p) - 1)],
        show_progress = show_progress,
        monodromy_options...,
    )
    additional_sols = solutions(additional_mres)
    if show_progress
        @info additional_mres
        @info "Computed $(length(additional_sols)) additional witnesses"
        @info "Compute trace using two parameter homotopies..."
    end

    # Perform parameter homotopy for the trace
    S = [map(s -> [s; 0], sols); additional_sols]
    γ = randn(ComplexF64)
    res1 = solve(
        verify_system,
        S;
        start_parameters = [0.0; base_params],
        target_parameters = [0.5 * γ; base_params],
        show_progress = show_progress,
        parameter_homotopy_options...,
    )
    S1 = solutions(res1)
    if length(S1) ≠ length(S)
        if show_progress
            @warn "Lost solution during parameter homotopy. Abort."
        end
        return nothing
    end

    res2 = solve(
        verify_system,
        solutions(res1);
        start_parameters = [0.5 * γ; base_params],
        target_parameters = [1.0 * γ; base_params],
        show_progress = show_progress,
        parameter_homotopy_options...,
    )
    S2 = solutions(res2)
    if length(S2) ≠ length(S)
        if show_progress
            @warn "Lost solution during parameter homotopy. Abort."
        end
        return nothing
    end

    # Compute trace now
    T = sum(S)
    T1 = sum(S1)
    T2 = sum(S2)

    M = [T T1 T2; 1 1 1]
    singvals = LA.svdvals(M)
    trace_norm = singvals[3] / singvals[1]

    if show_progress
        @info "Norm of trace: $trace_norm"
    end

    trace_norm < trace_tol
end
