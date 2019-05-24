export PathResult, PathTrackerStatus, PathTracker,
       pathtracker, pathtracker_startsolutions, solution,
       accuracy, residual, start_solution, isfailed, isatinfinity,
       issingular, isnonsingular, isprojective, isaffine, set_parameters!


const pathtracker_supported_keywords = [
    :at_infinity_check, :min_step_size_endgame_start,
    :min_val_accuracy, :samples_per_loop,
    :max_winding_number, :max_affine_norm,
    :overdetermined_min_accuracy, :overdetermined_min_residual]

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
        post_check_failed
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
    min_step_size_endgame_start::Float64
    min_val_accuracy::Float64
    samples_per_loop::Int
    max_winding_number::Int
    # Fallback tolerance only used if we track a path to the end without agreeing
    # on a valuation
    max_affine_norm::Float64
    # The minimal residual a solution needs to have to be considered
    # a solution of the original system (only applied for singular solutions)
    overdetermined_min_residual::Float64
    overdetermined_min_accuracy::Float64
end

function PathTrackerOptions(;
            at_infinity_check=true,
            min_step_size_endgame_start::Float64=1e-10,
            min_val_accuracy::Float64=0.001,
            samples_per_loop::Int=8,
            max_winding_number::Int=12,
            max_affine_norm::Float64=1e6,
            overdetermined_min_residual::Float64=1e-3,
            overdetermined_min_accuracy::Float64=1e-4)
    PathTrackerOptions(at_infinity_check, min_step_size_endgame_start, min_val_accuracy,
                       samples_per_loop, max_winding_number, max_affine_norm,
                       overdetermined_min_residual,
                       overdetermined_min_accuracy)
end

mutable struct PathTrackerState{V<:AbstractVector}
    status::PathTrackerStatus.states
    t::Float64
    prediction::V
    solution::V
    solution_accuracy::Float64
    solution_cond::Float64
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
    solution = copy(x)
    solution_accuracy = solution_cond = NaN
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

    PathTrackerState(status, t, prediction, solution, solution_accuracy, solution_cond,
                     endgame_zone_start, winding_number,
                     val, val_accuracy, prev_val, prev_val_accuracy)
end


function reset!(state::PathTrackerState)
    state.status = PathTrackerStatus.tracking
    state.t = 1.0
    state.prediction .= zero(eltype(state.prediction))
    state.solution .= zero(eltype(state.prediction))
    state.solution_accuracy = state.solution_cond = NaN
    state.endgame_zone_start = nothing
    state.winding_number = 0
    state.val .= 0.0
    state.val_accuracy .= Inf
    state.prev_val .= 0.0
    state.prev_val_accuracy .= Inf

    state
end


struct PathTrackerCache{T, V<:AbstractVector{Complex{T}}, S<:AbstractSystem, NC<:NewtonCache}
    unit_roots::Vector{ComplexF64}
    base_point::V
    target_system::S
    target_residual::Vector{Complex{T}}
    target_jacobian::Matrix{Complex{T}}
    target_newton_cache::NC
    weighted_ip::WeightedIP{T}
end

function PathTrackerCache(prob::Problem, core_tracker::CoreTracker)
    unit_roots = ComplexF64[]
    base_point = copy(core_tracker.state.x)

    x = currx(core_tracker)
    if is_squared_up_system(core_tracker.homotopy)
        # core_tracker.homotopy.target is a SquaredUpSystem
        F = core_tracker.homotopy.target.F
    else
        F = FixedHomotopy(core_tracker.homotopy, 0.0)
    end
    if affine_tracking(core_tracker)
        target_system = F
    elseif pull_back_is_to_affine(prob, x)
        target_system = PatchedSystem(F, state(EmbeddingPatch(), x))
    else # result is projective
        target_system = PatchedSystem(F, state(OrthogonalPatch(), x))
    end
    target_newton_cache = NewtonCache(target_system, x)

    res = evaluate(target_system, currx(core_tracker), target_newton_cache.system_cache)
    jac = jacobian(target_system, currx(core_tracker), target_newton_cache.system_cache)

    weighted_ip = WeightedIP(x)
    PathTrackerCache(unit_roots, base_point, target_system, res, jac, target_newton_cache, weighted_ip)
end

is_squared_up_system(::StraightLineHomotopy{<:Any,<:SquaredUpSystem}) = true
is_squared_up_system(::AbstractHomotopy) = false

"""
     PathTracker{Prob<:AbstractProblem, T, V<:AbstractVector{T}, CT<:CoreTracker}

`PathTracker` the way to track single paths. It combines the core path tracking routine
with an endgame, i.e., it can also deal with singular solutions as well as paths going to infinity.
`PathTracker` is a wrapper around [`CoreTracker`](@ref)
and thus has all configuration possibilities [`CoreTracker`](@ref) has.

There are the following `PathTracker` specific options (with their defaults in parens):

* `at_infinity_check::Bool=true`: Whether the path tracker should stop paths going to infinity early.
* `min_step_size_endgame_start=1e-10`: The endgame only starts if the step size becomes smaller that the provided value.
* `samples_per_loop::Int=5`: To compute singular solutions Cauchy's integral formula is used. The accuracy of the solutions increases with the number of samples per loop.
* `max_winding_number::Int=12`: The maximal number of loops used in Cauchy's integral formula.
* `max_affine_norm::Float64=1e6`: A fallback heuristic to decide whether a path is going to infinity.
* `min_val_accuracy::Float64=0.001`: A tolerance used to decide whether we are in the endgame zone.
* `overdetermined_min_accuracy=1e-5`: The minimal accuracy a non-singular solution needs to have to be considered a solution of the original system.
* `overdetermined_min_residual=1e-3`: The minimal residual a singular solution needs to have to be considered a solution of the original system.

In order to construct a pathtracker it is recommended to use the [`pathtracker`](@ref) and
[`pathtracker_startsolutions`](@ref) helper functions.
"""
struct PathTracker{V<:AbstractVector, Prob<:AbstractProblem, PTC<:PathTrackerCache, CT<:CoreTracker}
    problem::Prob
    core_tracker::CT
    state::PathTrackerState{V}
    options::PathTrackerOptions
    cache::PTC
end

function PathTracker(problem::AbstractProblem, core_tracker::CoreTracker; at_infinity_check=default_at_infinity_check(problem), optionskwargs...)
    state = PathTrackerState(core_tracker.state.x)
    options = PathTrackerOptions(; at_infinity_check=at_infinity_check, optionskwargs...)
    cache = PathTrackerCache(problem, core_tracker)
    PathTracker(problem, core_tracker, state, options, cache)
end

Base.show(io::IO, tracker::PathTracker) = print(io, "PathTracker")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathTracker) = x


default_at_infinity_check(prob::Problem{AffineTracking}) = true
default_at_infinity_check(prob::Problem{ProjectiveTracking}) = homvars(prob) !== nothing

"""
    track!(tracker::PathTracker, x₁, t₁::Float64=1.0; options...)::PathTrackerStatus.states

Track the path with start solution `x₁` from `t₁` towards `t=0`.

Possible values for the options are
* `accuracy::Float64`
* `max_corrector_iters::Int`
* `max_steps::Int`
* `start_parameters::AbstractVector`
* `target_parameters::AbstractVector`
"""
function track!(tracker::PathTracker, x₁, t₁::Float64=1.0;
        start_parameters::Union{Nothing,<:AbstractVector}=nothing,
        target_parameters::Union{Nothing,<:AbstractVector}=nothing,
        kwargs...)

    @unpack core_tracker, state, options, cache = tracker
    prev_options = set_options!(core_tracker; kwargs...)

    if start_parameters !== nothing
        set_start_parameters!(core_tracker.homotopy, start_parameters)
    end
    if target_parameters !== nothing
        set_target_parameters!(core_tracker.homotopy, target_parameters)
    end

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


        # Check at infinity
        if tracker.options.at_infinity_check &&
           core_tracker.state.Δs < options.min_step_size_endgame_start &&
           check_at_infinity(tracker)
            # We store when we entered the endgame zone
            if state.endgame_zone_start === nothing
                state.endgame_zone_start = t
            end

            state.status = PathTrackerStatus.at_infinity
            break
        end

        # Check that all valuations are accurate and solution is possibly singular
        if check_finite_singular_candidate(tracker; at_infinity_check=tracker.options.at_infinity_check)
            # We store when we entered the endgame zone
            if state.endgame_zone_start === nothing
                state.endgame_zone_start = t
            end

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
    if options.at_infinity_check && state.status == PathTrackerStatus.success
        # If the path is at infinity, then one of the homogenization variables is ≈ 0
        if (maximum(state.val_accuracy) < options.min_val_accuracy && minimum(state.val) < -0.05) ||
            vector_at_infinity(currx(core_tracker), options.max_affine_norm)
            state.status = PathTrackerStatus.at_infinity
        end
    end

    check_and_refine_solution!(tracker)
    # We have to set the number of blas threads to the previous value
    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)

    set_options!(core_tracker; prev_options...)
    state.status
end

"""
    track(tracker::PathTracker, x₁, t₁::Float64=1.0; path_number::Int=1, details::Symbol=:default, options...)::PathResult

Track the path with start solution `x₁` from `t₁` towards `t=0`. The `details` options controls
the level of details of the informations available in `PathResult`.

Possible values for the options are
* `accuracy::Float64`
* `max_corrector_iters::Int`
* `max_steps::Int`
* `start_parameters::AbstractVector`
* `target_parameters::AbstractVector`
"""
function track(tracker::PathTracker, x₁, t₁::Float64=1.0; path_number::Int=1, details::Symbol=:default, kwargs...)
    track!(tracker, x₁, t₁; kwargs...)
    PathResult(tracker, x₁, path_number; details=details)
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


"""
    set_parameters!(tracker::PathTracker; start_parameters=nothing, target_parameters=nothing)

Set the parameters of a parameter homotopy.
"""
function set_parameters!(tracker::PathTracker; start_parameters=nothing, target_parameters=nothing)
    if start_parameters !== nothing
        set_start_parameters!(tracker.core_tracker.homotopy, start_parameters)
    end
    if target_parameters !== nothing
        set_target_parameters!(tracker.core_tracker.homotopy, target_parameters)
    end
    nothing
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

function check_finite_singular_candidate(tracker::PathTracker; at_infinity_check=true)
    # We have a singular solution if one of the components
    # has a fractional valuation
    fractional_valuation = false
    for (i, ωᵢ) in enumerate(tracker.state.val)
        if !is_val_accurate(tracker.state, i, tracker.options) ||
            (at_infinity_check && ωᵢ ≤ -tracker.options.min_val_accuracy)
            return false
        end
        if abs(round(ωᵢ) - ωᵢ) > 0.1
            fractional_valuation = true
        end
    end

    fractional_valuation && return true

    # 3) We actually want some bad conditioning to not waste resources
    if unpack(tracker.core_tracker.state.jacobian.digits_lost, 0.0) > 4.0 ||
        tracker.core_tracker.state.ω > 100 ||
        tracker.core_tracker.state.Δs < 1e-4
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

"""
    check_and_refine_solution!(pathtracker)

In the case of success, we store the solution in `state.solution`
and try to refine the solution to the desired accuracy.
This also makes sure that our accuracy is correct after a possible pull back.
"""
function check_and_refine_solution!(tracker::PathTracker)
    @unpack core_tracker, state, options, cache = tracker

    if state.status ≠ PathTrackerStatus.success
        state.solution .= currx(core_tracker)
        return nothing
    end

    # state.winding_number > 0 only if during tracking we had the need to
    # use the cauchy endgame
    if state.winding_number > 0
        state.solution .= state.prediction
    else
        try
            res = refine!(core_tracker)
            state.solution .= currx(core_tracker)
        catch e
            # okay we had a singular solution after all
            if isa(e, LinearAlgebra.SingularException)
                tracker.state.winding_number = 1
                state.solution .= currx(core_tracker)
                for (ωᵢ, acc) in zip(state.prev_val, state.prev_val_accuracy)
                    if ωᵢ < -0.1 &&  acc < 0.05
                        state.status = PathTrackerStatus.at_infinity
                    end
                end
                if state.status != PathTrackerStatus.at_infinity &&
                   residual(tracker) > 0
                    state.status = PathTrackerStatus.at_infinity
                end
            else
                rethrow(e)
            end
        end
    end

    projective_refinement = false
    # Bring vector onto "standard form"
    if !affine_tracking(core_tracker)
        if pull_back_is_to_affine(tracker.problem, state.solution)
            ProjectiveVectors.affine_chart!(state.solution)
        else # projective result -> bring on (product of) sphere(s)
            LinearAlgebra.normalize!(state.solution)
            changepatch!(cache.target_system.patch, state.solution)
            projective_refinement = true
        end
    end

    # Catch the case that he tracker ran through but we still got a singular solution
    cond = condition_jacobian(tracker)
    if tracker.state.winding_number == 0 && cond > 1e13
        # make sure that we don't have a solution at infinity accidentaly
        for (ωᵢ, acc) in zip(state.prev_val, state.prev_val_accuracy)
            if ωᵢ < -0.1 &&  acc < 0.05
                state.status = PathTrackerStatus.at_infinity
            end
        end
        if state.status != PathTrackerStatus.at_infinity &&
           residual(tracker) > 0
            state.status = PathTrackerStatus.at_infinity
        end
        tracker.state.winding_number = 1
    elseif tracker.state.winding_number == 1 && cond < 1e5
        tracker.state.winding_number = 0
    end
    state.solution_cond = cond

    # If winding_number == 0 the path tracker simply ran through
    if (tracker.state.winding_number == 0)
        try
            # We need to reach the requested refinement_accuracy
            x̂ = cache.base_point
            if is_squared_up_system(core_tracker.homotopy)
                tol = options.overdetermined_min_accuracy
            else
                tol = core_tracker.options.refinement_accuracy
            end
            if projective_refinement
                result = newton!(x̂, cache.target_system, state.solution,
                                euclidean_norm, cache.target_newton_cache;
                                tol=tol, miniters=1, maxiters=2)
            else
                init_auto_scaling!(cache.weighted_ip, state.solution, AutoScalingOptions())
                result = newton!(x̂, cache.target_system, state.solution,
                                cache.weighted_ip, cache.target_newton_cache;
                                tol=tol, miniters=1, maxiters=2)
            end
            if isconverged(result)
                state.solution .= x̂
                if !pull_back_is_to_affine(tracker.problem, state.solution)
                    LinearAlgebra.normalize!(state.solution)
                    changepatch!(cache.target_system.patch, state.solution)
                end
                state.solution_accuracy = result.accuracy
            else
                state.solution_accuracy = result.accuracy
                state.status = PathTrackerStatus.post_check_failed
            end
        catch e
            # okay we had a singular solution after all
            if isa(e, LinearAlgebra.SingularException)
                tracker.state.winding_number = 1
            else
                rethrow(e)
            end
        end
    # For singular solutions check
    elseif is_squared_up_system(core_tracker.homotopy) &&
           residual(tracker) > options.overdetermined_min_residual
            state.status = PathTrackerStatus.post_check_failed
    end

    nothing
end

#############################
# Convencience constructors #
#############################

const pathtracker_startsolutions_supported_keywords = [
    problem_startsolutions_supported_keywords;
    coretracker_supported_keywords;
    pathtracker_supported_keywords]

"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracker`](@ref) and start solutions in the same way [`solve`](@ref) does it.
This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; system_scaling=true, kwargs...)
    invalid = invalid_kwargs(kwargs, pathtracker_startsolutions_supported_keywords)
    check_kwargs_empty(invalid, pathtracker_startsolutions_supported_keywords)

    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; system_scaling=system_scaling, supported...)
    core_tracker_supported, pathtracker_kwargs = splitkwargs(rest, coretracker_supported_keywords)
    core_tracker = CoreTracker(prob, start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); core_tracker_supported...)
    tracker = PathTracker(prob, core_tracker; pathtracker_kwargs...)
    (tracker=tracker, startsolutions=startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracker`](@ref) in the same way [`solve`](@ref) does it.
This also takes the same input arguments as `solve` with the exception that you do not need to specify startsolutions.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parameterized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = solution(track(tracker, x_a))
```
"""
pathtracker(args...; kwargs...) = first(pathtracker_startsolutions(args...; kwargs...))


############
## Result ##
############

solution(tracker::PathTracker) = pull_back(tracker.problem, tracker.state.solution)

function winding_number(tracker::PathTracker)
    tracker.state.winding_number == 0 ? nothing : tracker.state.winding_number
end

"""
    PathResult{V<:AbstractVector}

A `PathResult` is the result of tracking of a path using [`PathTracker`](@ref).
Its fields are

* `return_code`: One of `:success`, `:at_infinity` or any error code in [`PathTrackerStatus.states`](@ref) converted to a `Symbol`.
* `solution::V`: The solution vector.
* `t::Float64`: The value of `t` at which `solution` was computed. Note that if `return_code` is `:at_infinity`, then `t` is the value when this was decided.
* `accuracy::Union{Nothing, Float64}`: An approximation of ``||x-x^*||₂`` where ``x`` is the computed solution and ``x^*`` is the true solution.
* `residual::Union{Nothing, Float64}`: The value of the 2-norm of `H(solution, 0)`.
* `condition_jacobian::Union{Nothing, Float64}`: This is the condition number of the row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.
* `winding_number:Union{Nothing, Int}`: The estimated winding number. This is a lower bound on the multiplicity of the solution.
* `endgame_zone_start::Union{Nothing, Float64}`: The value of `t` at which we entered the endgame zone, i.e., where the path `x(t)` has an expansion a convergent Puiseux series near `t=0`.
* `start_solution::Union{Nothing, Int}`: The start solution of the path.
* `accepted_steps::Int`: The number of accepted steps during the path tracking.
* `rejected_steps::Int`: The number of rejected steps during the path tracking.
* `valuation::Union{Nothing, Vector{Float64}}`: An approximation of the valuation of the Puiseux series expansion of `x(t)`.
* `valuation_accuracy::Union{Nothing, Vector{Float64}}`: An estimate of the accuracy of the valuation of the Puiseux series expansion of `x(t)`.

     PathResult(tracker::PathTracker, start_solution=nothing, path_number::Union{Nothing,Int}=nothing; details=:default)

Possible `details` values are `:minimal` (minimal details), `:default` (default) and `:extensive` (all information possible).
"""
struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Union{Nothing, Float64}
    residual::Union{Nothing, Float64} # level 1+
    condition_jacobian::Union{Nothing, Float64}
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

function PathResult(tracker::PathTracker, start_solution, path_number::Union{Nothing,Int}=nothing; details::Symbol=:default)
    @unpack state, core_tracker, cache = tracker
    details_level = detailslevel(details)
    return_code = Symbol(state.status)
    windingnumber = winding_number(tracker)
    x = solution(tracker)
    # accuracy
    if !isnan(state.solution_accuracy)
        accuracy = state.solution_accuracy
    else
        accuracy = nothing
    end

    t = state.t
    # condition
    if isnan(state.solution_cond)
        condition_jac = nothing
    else
        condition_jac = state.solution_cond
    end
    # residual
    if return_code == :success && details_level ≥ 1
        res = residual(tracker)
    else
        res = nothing
    end

    endgame_zone_start = tracker.state.endgame_zone_start
    if details_level ≥ 1
        # mimic the behaviour in track! to get a start solution of the same type as x
        embed!(cache.base_point, tracker.problem, start_solution)
        startsolution = pull_back(tracker.problem, cache.base_point; regauge=false)
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
               windingnumber, endgame_zone_start, path_number,
               startsolution, accepted_steps, rejected_steps,
               valuation, valuation_accuracy)
end

function detailslevel(details::Symbol)
    if details == :minimal
        return 0
    elseif details == :extensive
        return 2
    else
        return 1
    end
end

function residual(tracker::PathTracker, x=tracker.state.solution)
    evaluate!(tracker.cache.target_residual, tracker.cache.target_system, x,
              tracker.cache.target_newton_cache.system_cache)
    euclidean_norm(tracker.cache.target_residual)
end

function condition_jacobian(tracker::PathTracker, x=tracker.state.solution)
    jacobian!(tracker.cache.target_jacobian, tracker.cache.target_system, x,
              tracker.cache.target_newton_cache.system_cache)
    row_scaling!(tracker.cache.target_jacobian)
    LinearAlgebra.cond(tracker.cache.target_jacobian)
end

"""
    result_type(tracker::PathTracker)

Returns the type of result `track` will return.
"""
result_type(tracker::PathTracker) = PathResult{typeof(solution(tracker))}

function Base.show(io::IO, r::PathResult)
    iscompact = get(io, :compact, false)
    if iscompact || haskey(io, :typeinfo)
        println(io, "• return_code: $(r.return_code)")
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

Get the accuracy of the solution ``x`` of the path, i.e., ``||H(x, 0)||₂``.
"""
accuracy(r::PathResult) = r.accuracy


"""
    residual(pathresult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, 0)||₂``.
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
    issingular(pathresult; tol=1e10)

Checks whether the path result is singular. This is true if
the winding number is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
issingular(r::PathResult; tol=1e10) = issingular(r, tol)
function issingular(r::PathResult, tol::Real)
    (unpack(r.winding_number, 0) ≥ 1 || unpack(r.condition_jacobian, 1.0) > tol) && LinearAlgebra.issuccess(r)
end

"""
    isnonsingular(pathresult; tol=1e10)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
isnonsingular(r::PathResult; kwargs...) = !issingular(r; kwargs...) && LinearAlgebra.issuccess(r)
isnonsingular(r::PathResult, tol::Real) = !issingular(r, tol) && LinearAlgebra.issuccess(r)


"""
    isreal(pathresult; tol=1e-6)

We consider a result as `real` if the 2-norm of the imaginary part of the solution is at most `tol`.
"""
Base.isreal(r::PathResult; tol=1e-6) = isreal(r, tol)
Base.isreal(r::PathResult, tol::Real) = isrealvector(r.solution, tol)

"""
    isprojective(pathresult)

Return`s true if the solution is a projective vector.
"""
isprojective(r::PathResult{<:PVector}) = true
isprojective(r::PathResult) = false

"""
    isaffine(pathresult)

Return`s true if the solution is an affine vector.
"""
isaffine(r::PathResult) = !isprojective(r)
