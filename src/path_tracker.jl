export PathResult, PathTrackerStatus, PathTracker,
       pathtracker, pathtracker_startsolutions, solution,
       accuracy, residual, start_solution, isfailed, isatinfinity,
       issingular, isnonsingular, isprojective


const pathtracker_supported_keywords = [
    :at_infinity_check, :max_step_size_endgame_start,
    :min_val_accuracy, :samples_per_loop,
    :max_winding_number, :max_affine_norm]

module PathTrackerStatus

    import ..CoreTrackerStatus

    @doc """
        PathTrackerStatus.states

    The possible return codes the path tracker can return are

    * `PathTrackerStatus.success`
    * `PathTrackerStatus.at_infinity`
    * `PathTrackerStatus.terminated_maximal_iterations`
    * `PathTrackerStatus.terminated_invalid_startvalue`
    * `PathTrackerStatus.terminated_step_size_too_small`
    * `PathTrackerStatus.terminated_singularity`
    * `PathTrackerStatus.terminated_ill_conditioned`
    * `PathTrackerStatus.terminated`
    """
    @enum states begin
        tracking
        success
        at_infinity
        terminated_invalid_startvalue
        terminated_maximal_iterations
        terminated_step_size_too_small
        terminated_singularity
        terminated_ill_conditioned
        tracker_failed
    end

    function status(code::CoreTrackerStatus.states)
        if code == CoreTrackerStatus.success
            return success
        elseif code == CoreTrackerStatus.terminated_invalid_startvalue
            return terminated_invalid_startvalue
        elseif code == CoreTrackerStatus.terminated_maximal_iterations
            return terminated_maximal_iterations
        elseif code == CoreTrackerStatus.terminated_step_size_too_small
            return terminated_step_size_too_small
        elseif code == CoreTrackerStatus.terminated_singularity
            return terminated_singularity
        elseif code == CoreTrackerStatus.terminated_ill_conditioned
            return terminated_ill_conditioned
        else
            # this shouldn't happen
            return tracker_failed
        end
    end
end

mutable struct PathTrackerOptions
    at_infinity_check::Bool
    max_step_size_endgame_start::Float64
    min_val_accuracy::Float64
    samples_per_loop::Int
    max_winding_number::Int
    # Fallback tolerance only used if we track a path to the end without agreeing
    # on a valuation
    max_affine_norm::Float64
end

function PathTrackerOptions(;
            at_infinity_check=true,
            max_step_size_endgame_start::Float64=1e-6,
            min_val_accuracy::Float64=0.001,
            samples_per_loop::Int=5,
            max_winding_number::Int=12,
            max_affine_norm::Float64=1e6)
    PathTrackerOptions(at_infinity_check, max_step_size_endgame_start, min_val_accuracy,
                       samples_per_loop, max_winding_number, max_affine_norm)
end

mutable struct PathTrackerState{V<:AbstractVector}
    status::PathTrackerStatus.states
    t::Float64
    prediction::V
    endgame_zone_start::Union{Nothing, Float64}
    winding_number::Int
    # Affine valuation
    val::Vector{Float64}
    val_accuracy::Vector{Float64}
    prev_val::Vector{Float64}
    prev_val_accuracy::Vector{Float64}
end

function PathTrackerState(x)
    status = PathTrackerStatus.tracking
    t = 1.0
    prediction = copy(x)
    endgame_zone_start = nothing
    winding_number = 0
    n = length(x)
    if x isa ProjectiveVectors.PVector
        n -= length(ProjectiveVectors.dims(x))
    end
    val = zeros(n)
    val_accuracy = copy(val)
    prev_val = copy(val)
    prev_val_accuracy = copy(val)

    PathTrackerState(status, t, prediction, endgame_zone_start, winding_number,
                     val, val_accuracy, prev_val, prev_val_accuracy)
end


function reset!(state::PathTrackerState)
    state.status = PathTrackerStatus.tracking
    state.t = 1.0
    state.prediction .= zero(eltype(state.prediction))
    state.endgame_zone_start = nothing
    state.winding_number = 0
    state.val .= 0.0
    state.val_accuracy .= Inf
    state.prev_val .= 0.0
    state.prev_val_accuracy .= Inf

    state
end


struct PathTrackerCache{T, V<:AbstractVector{T}, H<:HomotopyWithCache}
    unit_roots::Vector{ComplexF64}
    base_point::V
    residual::Vector{T}
    jacobian::Matrix{T}
    homotopy::H
end

function PathTrackerCache(core_tracker::CoreTracker)
    unit_roots = ComplexF64[]
    base_point = copy(core_tracker.state.x)

    x, t = currx(core_tracker), currt(core_tracker)
    if affine_tracking(core_tracker)
        homotopy = HomotopyWithCache(core_tracker.homotopy, x, t)
    else
        patched_homotopy = PatchedHomotopy(core_tracker.homotopy, state(EmbeddingPatch(), x))
        homotopy = HomotopyWithCache(patched_homotopy, x, t)
    end
    res = evaluate(homotopy, currx(core_tracker), currt(core_tracker))
    jac = jacobian(homotopy, currx(core_tracker), currt(core_tracker))

    PathTrackerCache(unit_roots, base_point, res, jac, homotopy)
end


"""
     PathTracker{Prob<:AbstractProblem, T, V<:AbstractVector{T}, CT<:CoreTracker}


The path tracker is *the* way to trace single paths. It combines the core path tracking routine
with an endgame. Thus, it decides whether a path is going to infinity, has a regular solution
or a singular solution.

It also handles any input transformations (e.g. affine -> projective -> affine)
such that the user gets out what he puts in.
"""
struct PathTracker{Prob<:AbstractProblem, T, V<:AbstractVector{T}, CT<:CoreTracker}
    problem::Prob
    core_tracker::CT
    state::PathTrackerState{V}
    options::PathTrackerOptions
    cache::PathTrackerCache{T, V}
end

function PathTracker(problem::AbstractProblem, core_tracker::CoreTracker; at_infinity_check=homvars(problem) !== nothing, optionskwargs...)
    state = PathTrackerState(core_tracker.state.x)
    options = PathTrackerOptions(; at_infinity_check=at_infinity_check, optionskwargs...)
    cache = PathTrackerCache(core_tracker)
    PathTracker(problem, core_tracker, state, options, cache)
end


function track!(tracker::PathTracker, x₁, t₁::Float64=1.0; kwargs...)
    @unpack core_tracker, state, options, cache = tracker
    prev_options = set_options!(core_tracker; kwargs...)

    # For performance reasons we single thread blas
    n_blas_threads = single_thread_blas()

    embed!(core_tracker.state.x, tracker.problem, x₁)
    setup!(core_tracker, core_tracker.state.x, t₁, 0.0)
    reset!(state)

    while state.status == PathTrackerStatus.tracking
        step!(core_tracker)
        t = state.t = real(currt(core_tracker))
        check_terminated!(core_tracker)

        if core_tracker.state.status ≠ CoreTrackerStatus.tracking
            state.status = PathTrackerStatus.status(core_tracker.state.status)
            break
        end

        # We only care if we moved forward
        core_tracker.state.last_step_failed && continue

        # Our endgame strategy is split in different stages
        # 1) During the path tracking we update an approximation
        #    of the valuation of x(t) since this has for t ≈ 0
        #    an expansion as a Puiseux series.
        #    The problem is that we do not a-priori know what t ≈ 0
        #    means. Thus, we have to find a criterion to check when
        #    an approximation is accurate.
        #    For this we require that:
        #         |val(xᵢ(t_k)) - val(xᵢ(t_{k-1}))| / (log(t_k) - log(t_{k-1}))
        #    is smaller than a given tolerance (`min_val_accuracy`)
        #
        # 2) If one of the valuation coordinates is sufficiently accurate
        #    we can decide whether this coordinate goes to infinity or not
        #    This happens in `check_at_infinity`
        #
        # 3) Now assume that val(xᵢ(t_k)) ≥ 0 for all i. Then we know
        #    that this converges to a finite solution. Now we have to distinguish
        #    between singular and regular solutions.
        #    We can handle singular by the Cauchy Endgame whereas for regular
        #    solutions we simply can track towards 0.

        # update valuation and associated data
        update_val!(state, core_tracker)

        # We want don't want to start the endgame too early
        core_tracker.state.Δs < options.max_step_size_endgame_start || continue

        # Check at infinity
        if tracker.options.at_infinity_check && check_at_infinity(tracker)
            state.status = PathTrackerStatus.at_infinity
            break
        end

        # Check that all valuations are accurate
        check_valuation_accurate(tracker) || continue

        # @show state.val state.val_accuracy
        # We store when we entered the endgame zone
        if state.endgame_zone_start === nothing
            state.endgame_zone_start = t
        end
        # Check singular
        if check_singular_canidate(tracker)
            retcode = predict_with_cauchy_integral_method!(state, core_tracker, options, cache)

            if retcode == :success
                state.status = PathTrackerStatus.success
            elseif retcode == :max_winding_number
                continue
            else # path tracker failed during the loops -> break
                state.status = PathTrackerStatus.tracker_failed
            end
        end

    end

    # It can happen that we jump too fast to a non-singular solution at infinity
    # For example, this happens with the Griewank & Osborne system (seed 130793) for path 3.
    # This only happens if we couldn't agrre on a valuation
    if state.status == PathTrackerStatus.success
        # If the path is at infinity, then one of the homogenization variables is ≈ 0
        if vector_at_infinity(currx(core_tracker), options.max_affine_norm)
            state.status = PathTrackerStatus.at_infinity
        end
    end

    if state.status == PathTrackerStatus.success && state.winding_number ≤ 1
        refine!(core_tracker)
    end

    # We have to set the number of blas threads to the previous value
    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)

    set_options!(core_tracker; prev_options...)
    state.status
end

function track(tracker::PathTracker, x₁, t₁::Float64=1.0; path_number::Int=1, details_level::Int=1, kwargs...)
    track!(tracker, x₁, t₁; kwargs...)
    PathResult(tracker, x₁, path_number; details_level=details_level)
end


function set_options!(core_tracker::CoreTracker;
                      accuracy::Float64=core_tracker.options.accuracy,
                      max_corrector_iters::Int=core_tracker.options.max_corrector_iters,
                      max_steps::Int=core_tracker.options.max_steps)
    core_tracker.options.accuracy = accuracy
    core_tracker.options.max_corrector_iters = max_corrector_iters
    core_tracker.options.max_steps = max_steps
    (accuracy=accuracy, max_corrector_iters=max_corrector_iters, max_steps=max_steps)
end



################
## VALUATIONS ##
################
"""
    update_val!(state::PathTrackerState, core_tracker::CoreTracker)

Update the valuation estimate of ``x(t)`` and assess it's accuracy.
"""
function update_val!(state::PathTrackerState, core_tracker::CoreTracker)
    @unpack x, ẋ = core_tracker.state
    t_k = abs(currt(core_tracker))
    compute_val!(state, x, ẋ, t_k)

    # Compute accuracy
    Δt = core_tracker.state.Δs_prev
    log_Δt_k = log1p(Δt / t_k)
    for i in eachindex(state.val)
        state.prev_val_accuracy[i] = state.val_accuracy[i]
        state.val_accuracy[i] = abs((state.prev_val[i] - state.val[i]) / log_Δt_k)
    end
    state
end
function compute_val!(state::PathTrackerState, x::AbstractVector, ẋ, t)
    for i in eachindex(state.val)
        state.prev_val[i] = state.val[i]
        state.val[i] = val(x, ẋ, t, i)
    end
    state
end
function compute_val!(state::PathTrackerState, x::PVector, ẋ, t)
    j = 1
    for (rᵢ, homᵢ) in ProjectiveVectors.dimension_indices_homvars(x)
        val_homᵢ = val(x, ẋ, t, homᵢ)
        for k in rᵢ
            state.prev_val[j] = state.val[j]
            state.val[j] = val(x, ẋ, t, k) - val_homᵢ
            j += 1
        end
    end
    state
end

val(x, ẋ, t, k) = t * re_dot(x[k], ẋ[k]) / abs2(x[k])

"Computes `real(z * conj(ż))`."
re_dot(z, ż) = begin x, y = reim(z); ẋ, ẏ = reim(ż); x * ẋ + y * ẏ end


############
## CHECKS ##
############

function check_valuation_accurate(tracker::PathTracker)
    for (i, ωᵢ) in enumerate(tracker.state.val)
        is_val_accurate(tracker.state, i, tracker.options) || return false
    end
    true
end

"""
    is_val_accurate(state, i, options)

Checks whether ``val(xᵢ(t))`` is accurate (as specified by the tolerance in `options`).
"""
function is_val_accurate(state::PathTrackerState, i::Int, options::PathTrackerOptions)
    # 1) require that the previous valuation accuracy is already enough
    state.prev_val_accuracy[i] < options.min_val_accuracy &&
    # We also want that the accuracy increase
    # or it should be much much more accurate
    (state.val_accuracy[i] < state.prev_val_accuracy[i]) ||
     state.val_accuracy[i] < options.min_val_accuracy^2
end

function check_at_infinity(tracker::PathTracker)
    for (i, ωᵢ) in enumerate(tracker.state.val)
        # Case 1 ωᵢ < 0 => |xⱼ(t)| -> ∞ for t -> 0
        if is_val_accurate(tracker.state, i, tracker.options) && ωᵢ < -0.05
            return true
        end
    end

    false
end

function check_singular_canidate(tracker::PathTracker)
    # this assumes that valuation is accurate and all valuations are > 0

    # We have a singular solution if one of the components
    # has a fractional valuation
    for (i, ωᵢ) in enumerate(tracker.state.val)
        if abs(round(ωᵢ) - ωᵢ) > 0.1 # This should be sufficient for the most common cases
            return true
        end
    end

    # 3) We actually want some bad conditioning to not waste resources
    if tracker.core_tracker.state.digits_lost > 4 ||
        tracker.core_tracker.state.ω > 100 ||
        tracker.core_tracker.state.Δs < 1e-6
        return true
    end

    false
end

vector_at_infinity(x::AbstractVector, tol) = false
vector_at_infinity(x::PVector, tol) = ProjectiveVectors.norm_affine_chart(x) > tol


##########
# Cauchy #
##########

"""
    predict_with_cauchy_integral_method!(state, core_tracker, options, cache)

Try to predict the value of `x(0)` using [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula)
At each iteration we are at some point ``(x, t)``. We then track the polygon defined
by ``te^{i2πk/n}`` until we end again at ``x``. Here ``n`` is the number of samples we take
per loop.
Returns a symbol indicating whether the prediction was successfull.
"""
function predict_with_cauchy_integral_method!(state, core_tracker, options, cache)
    @unpack prediction = state
    @unpack unit_roots, base_point = cache
    @unpack samples_per_loop, max_winding_number = options

    compute_unit_roots!(unit_roots, samples_per_loop)

    t = real(currt(core_tracker))

    base_point .= currx(core_tracker)
    prediction .= zero(eltype(state.prediction))

    # during the loop we fix the affine patch
    @unpack accepted_steps, rejected_steps = core_tracker.state
    fix_patch!(core_tracker)

    m = k = 1
    while m ≤ max_winding_number
        θⱼ = t
        for j=1:samples_per_loop
            θⱼ₋₁ = θⱼ
            θⱼ = t * unit_roots[j]

            retcode = track!(core_tracker, currx(core_tracker), θⱼ₋₁, θⱼ)

            accepted_steps += core_tracker.state.accepted_steps
            rejected_steps += core_tracker.state.rejected_steps

            if retcode != CoreTrackerStatus.success
                # during the loop we fixed the affine patch
                unfix_patch!(core_tracker)
                @pack! core_tracker.state = accepted_steps, rejected_steps
                return Symbol(retcode)
            end

            prediction .+= currx(core_tracker)
        end

        if euclidean_distance(base_point, currx(core_tracker)) < 4 * core_tracker.options.accuracy
            break
        end

        m += 1
    end

    # we have to undo the fixing of the patch
    unfix_patch!(core_tracker)
    @pack! core_tracker.state = accepted_steps, rejected_steps

    if m > max_winding_number
        return :max_winding_number
    end

    state.winding_number = m
    prediction ./= m * samples_per_loop
    :success
end

fix_patch!(tracker::CoreTracker) = tracker.options.update_patch = false
unfix_patch!(tracker::CoreTracker) = tracker.options.update_patch = true

"""
    compute_unit_roots!(unit_roots, n)

Fill `unit_roots` with the unit roots ``exp(2πi/n)`` for ``i=1,...,n``.
"""
function compute_unit_roots!(unit_roots, n)
    if length(unit_roots) ≠ n
        resize!(unit_roots, n)
        unit_roots[1] = primitive_root = cis(2π/n)
        unit_roots[n] = complex(1.0)
        for k=2:(n-1)
            unit_roots[k] = unit_roots[k-1] * primitive_root
        end
    end
    unit_roots
end


"""
    type_of_x(pathtracker)

Returns the type of `x`.
"""
type_of_x(tracker::PathTracker) = typeof(tracker.core_tracker.state.x)

#############################
# Convencience constructors #
#############################

const pathtracker_startsolutions_supported_keywords = [
    problem_startsolutions_supported_keywords;
    coretracker_supported_keywords;
    pathtracker_supported_keywords]

"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracker`](@ref) and `startsolutions` in the same way `solve` does it.
This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    invalid = invalid_kwargs(kwargs, pathtracker_startsolutions_supported_keywords)
    check_kwargs_empty(invalid, pathtracker_startsolutions_supported_keywords)

    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    core_tracker_supported, pathtracker_kwargs = splitkwargs(rest, coretracker_supported_keywords)
    core_tracker = CoreTracker(prob, start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); rest...)
    tracker = PathTracker(prob, core_tracker; pathtracker_kwargs...)
    (tracker=tracker, startsolutions=startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracker`](@ref) in the same way `solve` does it.
This also takes the same input arguments as `solve` with the exception that you do not need to specify startsolutions.
This is convenient if you want to investigate single paths.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parameterized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = track(tracker, x_a).x
```

### Trace a path
To trace a path you can use the [`iterator`](@ref) method.

```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```
"""
function pathtracker(args...; kwargs...)
    tracker, _ = pathtracker_startsolutions(args...; kwargs...)
    tracker
end


############
## Result ##
############

function solution(tracker::PathTracker)
    if tracker.state.winding_number == 0
        pull_back(tracker.problem, currx(tracker.core_tracker))
    else
        pull_back(tracker.problem, tracker.state.prediction)
    end
end

function winding_number(tracker::PathTracker)
    tracker.state.winding_number == 0 ? nothing : tracker.state.winding_number
end

"""
     PathResult(tracker::PathTracker, start_solution=nothing, path_number::Union{Nothing,Int}=nothing; details_level=1)

Possible `details_levels` values are `0` (minimal details), `1` (default) and `2` (all information possible).
"""
struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Union{Nothing, Float64}
    residual::Union{Nothing, Float64} # level 1+
    condition_jacobian::Union{Nothing, Float64} # level 1+
    winding_number::Union{Nothing, Int}
    endgame_zone_start::Union{Nothing, Float64}
    path_number::Union{Nothing, Int}
    start_solution::Union{Nothing, V} # level 1+
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    valuation::Union{Nothing, Vector{Float64}} # level 2+
    valuation_accuracy::Union{Nothing, Vector{Float64}} # level 2+
end

function PathResult(tracker::PathTracker, start_solution, path_number::Union{Nothing,Int}=nothing; details_level::Int=1)
    @unpack state, core_tracker, cache = tracker

    return_code = Symbol(state.status)
    windingnumber = winding_number(tracker)
    x = solution(tracker)
    # accuracy
    if windingnumber === nothing && return_code == :success
        accuracy = core_tracker.state.accuracy
    else
        accuracy = nothing
    end

    t = state.t
    # residual
    if return_code == :success && details_level ≥ 1
        if affine_tracking(tracker.core_tracker)
            res = residual(tracker, currx(core_tracker), t)
            condition_jac = condition_jacobian(tracker, currx(core_Tracker), t)
        else
            # we simulate the affine vector
            # We cannot use `x` from above since pull_back also takes care of possible
            # permutations of the variables.
            cache.base_point .= currx(core_tracker)
            y = ProjectiveVectors.affine_chart!(cache.base_point)

            res = residual(tracker, y, t)
            condition_jac = condition_jacobian(tracker, y, t)
        end
    else
        res = nothing
        condition_jac = nothing
    end

    endgame_zone_start = tracker.state.endgame_zone_start
    if details_level ≥ 1
        # mimic the behaviour in track! to get a start solution of the same type as x
        startsolution = pull_back(tracker.problem, embed!(cache.base_point, tracker.problem, start_solution))
    else
        startsolution = nothing
    end

    accepted_steps = core_tracker.state.accepted_steps
    rejected_steps = core_tracker.state.rejected_steps

    if details_level == 2
        valuation = copy(tracker.state.val)
        valuation_accuracy = copy(tracker.state.val_accuracy)
    else
        valuation = nothing
        valuation_accuracy = nothing
    end

    PathResult(return_code, x, t, accuracy, res, condition_jac,
                      windingnumber, endgame_zone_start,
                      path_number,
                      startsolution,
                      accepted_steps, rejected_steps,
                      valuation, valuation_accuracy)
end

function residual(tracker::PathTracker, x, t)
    evaluate!(tracker.cache.residual, tracker.cache.homotopy, x, t)
    euclidean_norm(tracker.cache.residual)
end

function condition_jacobian(tracker, x, t)
    jacobian!(tracker.cache.jacobian, tracker.cache.homotopy, x, t)
    row_scaling!(tracker.cache.jacobian)
    LinearAlgebra.cond(tracker.cache.jacobian)
end



"""
    result_type(tracker::PathTracker)

Returns the type of result `track` will return.
"""
result_type(tracker::PathTracker) = PathResult{typeof(solution(tracker))}

function Base.show(io::IO, r::PathResult)
    iscompact = get(io, :compact, false)
    if iscompact || haskey(io, :typeinfo)
        println(io, " • return_code: $(r.return_code)")
        println(io, " • solution: ", r.solution)
        r.accuracy !== nothing &&
            println(io, " • accuracy: $(Printf.@sprintf "%.3e" r.accuracy)")
        r.winding_number !== nothing &&
            println(io, " • winding_number: $(r.winding_number)")
        r.path_number !== nothing &&
            println(io, " • path_number: ", r.path_number)
    else
        println(io, "PathResult")
        println(io, "=================")
        println(io, " • return_code: $(r.return_code)")
        println(io, " • solution: ", r.solution)
        r.accuracy !== nothing &&
            println(io, " • accuracy: $(Printf.@sprintf "%.3e" r.accuracy)")
        r.residual !== nothing &&
            println(io, " • residual: $(Printf.@sprintf "%.3e" r.residual)")
        r.winding_number !== nothing &&
            println(io, " • winding_number: $(r.winding_number)")
        r.condition_jacobian !== nothing &&
            println(io, " • condition_jacobian: $(Printf.@sprintf "%.3e" r.condition_jacobian)")
        r.path_number !== nothing &&
            println(io, " • path_number: ", r.path_number)
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathResult) = x

"""
    solution(pathresult)

Get the solution of the path.
"""
solution(r::PathResult) = r.solution


"""
    accuracy(pathresult)

Get the accuracy of the solution ``x`` of the path, i.e., ``||H(x, 0)||_2``.
"""
accuracy(r::PathResult) = r.accuracy


"""
    residual(pathresult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, 0)||_2``.
"""
residual(r::PathResult) = r.residual

"""
    start_solution(pathresult)

Get the start solution of the solution ``x`` of the path.
"""
start_solution(r::PathResult) = r.start_solution

"""
    issuccess(pathresult)

Checks whether the path is successfull.
"""
LinearAlgebra.issuccess(r::PathResult) = r.return_code == :success

"""
    isfailed(pathresult)

Checks whether the path failed.
"""
isfailed(r::PathResult) =!(r.return_code == :at_infinity || r.return_code == :success)


"""
    isatinfinity(pathresult)

Checks whether the path goes to infinity.
"""
isatinfinity(r::PathResult) = r.return_code == :at_infinity

"""
    isfinite(pathresult)

Checks whether the path result is finite.
"""
Base.isfinite(r::PathResult) = r.return_code == :success # we don't check isaffine to make other code easier

"""
    issingular(pathresult; tol=1e14)

Checks whether the path result is singular. This is true if
the winding number is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
issingular(r::PathResult; tol=1e14) = issingular(r, tol)
function issingular(r::PathResult, tol::Real)
    (unpack(r.winding_number, 0) > 1 || unpack(r.condition_jacobian, 1.0) > tol) && LinearAlgebra.issuccess(r)
end

"""
    isnonsingular(pathresult; tol=1e14)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
isnonsingular(r::PathResult; tol=1e14) = isnonsingular(r, tol)
isnonsingular(r::PathResult, tol::Real) = !issingular(r, tol) && LinearAlgebra.issuccess(r)


"""
    isreal(pathresult; tol=1e-6)

We consider a result as `real` if the 2-norm of the imaginary part of the solution is at most `tol`.
"""
Base.isreal(r::PathResult; tol=1e-6) = isreal(r, tol)
Base.isreal(r::PathResult, tol::Real) = isrealvector(r.solution, tol)

isprojective(r::PathResult{<:PVector}) = true
isprojective(r::PathResult) = false
