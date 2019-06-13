export PathResult, PathTrackerStatus, PathTracker,
       pathtracker, pathtracker_startsolutions, solution,
       accuracy, residual, start_solution, isfailed, isatinfinity,
       issingular, isnonsingular, isprojective, isaffine, set_parameters!


const pathtracker_supported_keywords = [
    :at_infinity_check, :min_step_size_endgame_start,
    :min_val_accuracy, :samples_per_loop,
    :max_winding_number, :max_affine_norm,
    :overdetermined_min_accuracy, :overdetermined_min_residual]


###############
## VALUATION ##
###############

struct Valuation
    v::Vector{Float64}
    v̇::Vector{Float64}
    v̈::Vector{Float64}
    #  0.5|x(s)|^2
    abs_x_sq::Vector{Float64}
    # d/ds 0.5|x(s)|^2
    abs_x_sq_dot::Vector{Float64}
    # d^2/ds^2 0.5|x(s)|^2
    abs_x_sq_ddot::Vector{Float64}
    # d^3/ds^3 0.5|x(s)|^2
    abs_x_sq_dddot::Vector{Float64}
end

@enum ValuationVerdict begin
    VALUATION_INDECISIVE
    VALUATION_FINITE
    VALUATION_AT_INFINIY
end

function Valuation(x::ProjectiveVectors.PVector; at_infinity_check::Bool=true)
    if at_infinity_check
        Valuation(length(x) - length(ProjectiveVectors.dims(x)))
    else
        Valuation(length(x))
    end
end
Valuation(x::Vector; kwargs...) = Valuation(length(x))
function Valuation(n::Integer)
    v = zeros(n)
    v̇ = zeros(n)
    v̈ = zeros(n)
    abs_x_sq = zeros(n)
    abs_x_sq_dot = zeros(n)
    abs_x_sq_ddot = zeros(n)
    abs_x_sq_dddot = zeros(n)
    Valuation(v, v̇, v̈, abs_x_sq, abs_x_sq_dot, abs_x_sq_ddot, abs_x_sq_dddot)
end

Base.show(io::IO, val::Valuation) = print_fieldnames(io, val)
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

function reset!(val::Valuation)
    val.v .= 0.0
    val.v̇ .= 0.0
    val.v̈ .= 0.0
    val.abs_x_sq .= 0.0
    val.abs_x_sq_dot .= 0.0
    val.abs_x_sq_ddot .= 0.0
    val.abs_x_sq_dddot .= 0.0
    val
end

function update!(val::Valuation, z::Vector, ż, Δs::Float64, ::Val{true})
    _update!(val, z, ż, Δs)
end
function update!(val::Valuation, z::PVector, ż, Δs::Float64, ::Val{false})
    _update!(val, z, ż, Δs)
end
function update!(val::Valuation, z::PVector, ż, Δs::Float64, ::Val{true})
    _update_affine!(val, z, ż, Δs)
end
function _update!(val::Valuation, z::AbstractVector, ż, Δs::Float64)
    v, v̇, v̈ = val.v, val.v̇, val.v̈
    for i in eachindex(z)
        x, y = reim(z[i])
        ẋ, ẏ = reim(ż[i])
        x²_plus_y² = abs2(z[i])
        vᵢ, v̇ᵢ = v[i], v̇[i]
        v[i] = -(x*ẋ + y*ẏ) / x²_plus_y²
        v̇[i] = (v[i] - vᵢ) / Δs
        v̈[i] = (v̇[i] - v̇ᵢ) / Δs

        abs_x_sq_dotᵢ, abs_x_sq_ddotᵢ = val.abs_x_sq_dot[i], val.abs_x_sq_ddot[i]

        val.abs_x_sq[i] = 0.5x²_plus_y²
        val.abs_x_sq_dot[i] = x*ẋ + y*ẏ
        val.abs_x_sq_ddot[i] = (val.abs_x_sq_dot[i] - abs_x_sq_dotᵢ) / Δs
        val.abs_x_sq_dddot[i] = (val.abs_x_sq_ddot[i] - abs_x_sq_ddotᵢ) / Δs
    end

    val
end

function _update_affine!(val::Valuation, z::PVector, ż, Δs::Float64)
    i = 1
    v, v̇, v̈ = val.v, val.v̇, val.v̈
    for (rⱼ, homⱼ) in ProjectiveVectors.dimension_indices_homvars(z)
        xⱼ, yⱼ = reim(z[homⱼ])
        ẋⱼ, ẏⱼ = reim(ż[homⱼ])
        xⱼ²_plus_yⱼ² = xⱼ^2 + xⱼ^2
        vⱼ = -(xⱼ*ẋⱼ + yⱼ*ẏⱼ) / xⱼ²_plus_yⱼ²
        for k in rⱼ
            x, y = reim(z[k])
            ẋ, ẏ = reim(ż[k])
            x²_plus_y² = x^2 + y^2
            vᵢ, v̇ᵢ = v[i], v̇[i]
            v[i] = -(x*ẋ + y*ẏ) / x²_plus_y²

            v[i] -= vⱼ
            v̇[i] = (v[i] - vᵢ) / Δs
            v̈[i] = (v̇[i] - v̇ᵢ) / Δs

            abs_x_sqᵢ, abs_x_sq_dotᵢ, abs_x_sq_ddotᵢ =
                                val.abs_x_sq[i], val.abs_x_sq_dot[i], val.abs_x_sq_ddot[i]
            val.abs_x_sq[i] = 0.5 * x²_plus_y² / xⱼ²_plus_yⱼ²
            val.abs_x_sq_dot[i] =  (val.abs_x_sq[i] - abs_x_sqᵢ) / Δs
            val.abs_x_sq_ddot[i] = (val.abs_x_sq_dot[i] - abs_x_sq_dotᵢ) / Δs
            val.abs_x_sq_dddot[i] = (val.abs_x_sq_ddot[i] - abs_x_sq_ddotᵢ) / Δs
            i += 1
        end
    end
    val
end

"""
    judge(val::Valuation; s_max::Float64=-log(eps()))::ValuationVerdict

Judge the current valuation.
"""
function judge(val::Valuation, J::Jacobian, s, s_max::Float64=-log(eps()))
    @unpack v, v̇, v̈, abs_x_sq, abs_x_sq_dot, abs_x_sq_ddot, abs_x_sq_dddot = val
    Δs = abs(s_max - s)
    status = VALUATION_INDECISIVE
    all_positive = true
    indecisive = false
    for i in eachindex(v)
        Δv̂ᵢ = Δs * abs(v̇[i]) + Δs^2 * abs(v̈[i])
        if abs(v̇[i]) > 1e-2
            indecisive = true
        end

        if status != VALUATION_AT_INFINIY &&
            v[i] < -1/20 &&
            Δv̂ᵢ < 1e-5 &&
            abs_x_sq[i] > -2v[i] * abs_x_sq_dot[i] > 0 &&
            abs_x_sq_dot[i] > -2v[i] * abs_x_sq_ddot[i] > 0.1 &&
            abs_x_sq_ddot[i] > -2v[i] * abs_x_sq_dddot[i] > 0.1

            status = VALUATION_AT_INFINIY
        end

        if v[i] < -0.05
            all_positive = false
        end
    end
    status == VALUATION_AT_INFINIY && return VALUATION_AT_INFINIY
    indecisive && return VALUATION_INDECISIVE
    # require s > 0 otherwise we can get false singular solutions
    all_positive && Δs > 0 && return VALUATION_FINITE

    VALUATION_INDECISIVE
end

@inline function update!(val::Valuation, core_tracker::CoreTracker; at_infinity_check::Bool=true)
    z = core_tracker.state.x
    ż = core_tracker.state.ẋ
    Δs = core_tracker.state.Δs_prev
    if at_infinity_check
        update!(val, z, ż, Δs, Val(true))
    else
        update!(val, z, ż, Δs, Val(false))
    end
end

function at_infinity_post_check(val::Valuation)
    for (vᵢ, v̇ᵢ) in zip(val.v, val.v̇)
        if vᵢ < -1/20 && abs(v̇ᵢ) < 1e-3
            return true
        end
    end
    false
end
#################
## PATHTRACKER ##
#################

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
    * `PathTrackerStatus.post_check_failed`
    * `PathTrackerStatus.excess_solution`
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
        excess_solution
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
            samples_per_loop::Int=12,
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
    s::Float64
    prediction::V
    solution::V
    solution_accuracy::Float64
    solution_cond::Float64
    endgame_zone_start::Union{Nothing, Float64}
    winding_number::Int
    val::Valuation
end

function PathTrackerState(x; at_infinity_check::Bool=true)
    status = PathTrackerStatus.tracking
    s = 0.0
    prediction = copy(x)
    solution = copy(x)
    solution_accuracy = solution_cond = NaN
    endgame_zone_start = nothing
    winding_number = 0
    val = Valuation(x; at_infinity_check=at_infinity_check)

    PathTrackerState(status, s, prediction, solution,
                    solution_accuracy, solution_cond,
                     endgame_zone_start, winding_number, val)
end


function reset!(state::PathTrackerState)
    state.status = PathTrackerStatus.tracking
    state.s = 0.0
    state.prediction .= zero(eltype(state.prediction))
    state.solution .= zero(eltype(state.prediction))
    state.solution_accuracy = state.solution_cond = NaN
    state.endgame_zone_start = nothing
    state.winding_number = 0
    reset!(state.val)

    state
end


struct PathTrackerCache{T, V<:AbstractVector{Complex{T}}, S<:AbstractSystem, NC<:AbstractNewtonCache}
    unit_roots::Vector{ComplexF64}
    base_point::V
    target_system::S
    target_residual::Vector{Complex{T}}
    target_jacobian::Jacobian{Complex{T}}
    target_newton_cache::NC
    weighted_ip::WeightedIP{T}
end

function PathTrackerCache(prob::Problem, core_tracker::CoreTracker)
    unit_roots = ComplexF64[]
    base_point = copy(core_tracker.state.x)

    x = currx(core_tracker)
    if is_squared_up_system(core_tracker.homotopy)
        F = squared_up_system(core_tracker.homotopy).F
    else
        F = FixedHomotopy(core_tracker.homotopy, Inf)
    end
    if affine_tracking(core_tracker)
        target_system = F
    elseif pull_back_is_to_affine(prob, x)
        target_system = PatchedSystem(F, state(EmbeddingPatch(), x))
    else # result is projective
        target_system = PatchedSystem(F, state(OrthogonalPatch(), x))
    end
    target_newton_cache = newton_cache(target_system, x)

    res = evaluate(target_system, currx(core_tracker), target_newton_cache.system_cache)
    jac = similar(res, size(target_system))
    weighted_ip = WeightedIP(x)
    PathTrackerCache(unit_roots, base_point, target_system, res, Jacobian(Random.rand!(jac)), target_newton_cache, weighted_ip)
end

is_squared_up_system(::StraightLineHomotopy{<:Any,<:SquaredUpSystem}) = true
is_squared_up_system(H::LogHomotopy) = is_squared_up_system(H.homotopy)
is_squared_up_system(::AbstractHomotopy) = false

squared_up_system(H::StraightLineHomotopy{<:Any,<:SquaredUpSystem}) = H.target
squared_up_system(H::LogHomotopy) = squared_up_system(H.homotopy)

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

function PathTracker(prob::AbstractProblem, x::AbstractVector{<:Number}; at_infinity_check=default_at_infinity_check(prob), min_step_size=eps(), kwargs...)
    core_tracker_supported, optionskwargs = splitkwargs(kwargs, coretracker_supported_keywords)
    core_tracker = CoreTracker(prob, x, complex(0.0), 36.0;
                        log_transform=true, predictor=Pade21(),
                        min_step_size=min_step_size,
                        core_tracker_supported...)
    state = PathTrackerState(core_tracker.state.x; at_infinity_check=at_infinity_check)
    options = PathTrackerOptions(; at_infinity_check=at_infinity_check, optionskwargs...)
    cache = PathTrackerCache(prob, core_tracker)
    PathTracker(prob, core_tracker, state, options, cache)
end

Base.show(io::IO, tracker::PathTracker) = print(io, "PathTracker")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathTracker) = x


seed(PT::PathTracker) = PT.problem.seed
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
function track!(tracker::PathTracker, x₁, s₁=0.0, s₀=36.0;
        start_parameters::Union{Nothing,<:AbstractVector}=nothing,
        target_parameters::Union{Nothing,<:AbstractVector}=nothing,
        kwargs...)

    @unpack core_tracker, state, options, cache = tracker
    prev_options = set_options!(core_tracker; kwargs...)
    set_parameters!(tracker; start_parameters=start_parameters, target_parameters=target_parameters)

    _track!(tracker, x₁, s₁, s₀)

    set_options!(core_tracker; prev_options...)
    state.status
end

function _track!(tracker::PathTracker, x₁, s₁::Real, s₀::Real)
    @unpack core_tracker, state, options, cache = tracker
    # For performance reasons we single thread blas
    n_blas_threads = single_thread_blas()

    embed!(core_tracker.state.x, tracker.problem, x₁)
    setup!(core_tracker, core_tracker.state.x, s₁, s₀)
    reset!(state)

    while state.status == PathTrackerStatus.tracking
        step!(core_tracker)
        state.s = real(currt(core_tracker))
        check_terminated!(core_tracker)
        if core_tracker.state.status != CoreTrackerStatus.tracking &&
           core_tracker.state.status != CoreTrackerStatus.success
            state.status = PathTrackerStatus.status(core_tracker.state.status)
            break
        end

        if core_tracker.state.jacobian.corank > 0
            state.status = PathTrackerStatus.terminated_ill_conditioned
            break
        end

        # We only care if we moved forward
        core_tracker.state.last_step_failed && continue

        # Our endgame strategy is split in different stages
        # 1) During the path tracking we update an approximation
        #    of the valuation of x(t) since this has for t ≈ 0
        #    an expansion as a Puiseux series.
        #
        # 2) If one of the valuation coordinates is sufficiently accurate
        #    we can decide whether this coordinate goes to infinity or not
        #    This happens in `judge`.
        #
        # 3) Now assume that val(xᵢ(t_k)) ≥ 0 for all i. Then we know
        #    that this converges to a finite solution. Now we have to distinguish
        #    between singular and regular solutions.
        #    We can handle singular by the Cauchy Endgame whereas for regular
        #    solutions we simply can track towards ∞.

        # update valuation and associated data
        update!(state.val, core_tracker; at_infinity_check=tracker.options.at_infinity_check)

        verdict = judge(state.val, core_tracker.state.jacobian, state.s, max(36.0, s₀))
        if tracker.options.at_infinity_check && verdict == VALUATION_AT_INFINIY
            state.status = PathTrackerStatus.at_infinity
            break
        end

        if verdict == VALUATION_FINITE && is_singularish(state.val, core_tracker)
            retcode = predict_with_cauchy_integral_method!(state, core_tracker, options, cache)
            if retcode == :success
                state.status = PathTrackerStatus.success
            elseif retcode == :max_winding_number
                continue
            else # path tracker failed during the loops -> break
                state.status = PathTrackerStatus.tracker_failed
            end
        end

        if core_tracker.state.status == CoreTrackerStatus.success
            state.status = PathTrackerStatus.success
            break
        end
    end

    check_and_refine_solution!(tracker)
    # We have to set the number of blas threads to the previous value
    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)
    state.status
end


function is_singularish(val::Valuation, core_tracker::CoreTracker)
    # 1) check if val is rational for some i, then the solution is singular
    for v in val.v
        v < -0.05 && return false # this function shouldn't be called in the first
        abs(round(v) - v) > 0.1 && return true
    end
    # 2) use condition number / digits_lost to decide
    unpack(core_tracker.state.jacobian.cond, 0.0) > 1e6 ||
    unpack(core_tracker.state.jacobian.digits_lost, 0.0) > 4.0
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
function track(tracker::PathTracker, x₁, t₁=0.0, t₀=36.0; path_number::Int=1, details::Symbol=:default, kwargs...)
    track!(tracker, x₁, t₁, t₀; kwargs...)
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
@inline function set_parameters!(tracker::PathTracker; start_parameters=nothing, target_parameters=nothing)
    if start_parameters !== nothing
        set_start_parameters!(basehomotopy(tracker.core_tracker.homotopy), start_parameters)
    end
    if target_parameters !== nothing
        set_target_parameters!(basehomotopy(tracker.core_tracker.homotopy), target_parameters)
    end
    nothing
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

    # compute_unit_roots!(unit_roots, samples_per_loop)
    initial_step_size = core_tracker.options.initial_step_size
    initial_max_steps = core_tracker.options.max_steps
    s = real(currt(core_tracker))

    base_point .= currx(core_tracker)
    prediction .= zero(eltype(state.prediction))

    # during the loop we fix the affine patch
    @unpack accepted_steps, rejected_steps = core_tracker.state


    fix_patch!(core_tracker)

    m = k = 1
    ∂θ = 2π / samples_per_loop
    core_tracker.options.initial_step_size = ∂θ #0.5∂θ
    core_tracker.options.max_steps = 25
    while m ≤ max_winding_number
        θⱼ = 0.0
        for j=1:samples_per_loop
            θⱼ₋₁ = θⱼ
            θⱼ += ∂θ

            retcode = track!(core_tracker, currx(core_tracker), s + im*θⱼ₋₁, s + im*θⱼ; checkstartvalue=false)
            accepted_steps += core_tracker.state.accepted_steps
            rejected_steps += core_tracker.state.rejected_steps

            if retcode != CoreTrackerStatus.success
                # during the loop we fixed the affine patch
                unfix_patch!(core_tracker)
                core_tracker.options.initial_step_size = initial_step_size
                core_tracker.options.max_steps = initial_max_steps
                @pack! core_tracker.state = accepted_steps, rejected_steps

                return Symbol(retcode)
            end

            prediction .+= currx(core_tracker)
        end

        if euclidean_distance(base_point, currx(core_tracker)) < 4core_tracker.options.accuracy
            break
        end

        m += 1
    end
    core_tracker.options.initial_step_size = initial_step_size
    core_tracker.options.max_steps = initial_max_steps
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

    # The tracking failed. Let's look at the last valuation and see whether we can
    # resonably classify this as going to infinity
    if (state.status == PathTrackerStatus.terminated_ill_conditioned ||
        state.status == PathTrackerStatus.terminated_step_size_too_small ||
        state.status == PathTrackerStatus.terminated_singularity) &&
        at_infinity_post_check(state.val)

       state.status = PathTrackerStatus.at_infinity
   end

    if state.status ≠ PathTrackerStatus.success
        state.solution .= currx(core_tracker)
        return nothing
    end

    # We have to differentiate

    # 1.a) Non-singular solution of regular system
    # 1.b) Singular solution of regular system (Cauchy engame used)
    # 2.a) Non-singular solution of squared up system
    # 2.b) Singular solution of squared up system (Cauchy engame used)

    is_singular = state.winding_number > 0
    if is_singular
        state.solution .= state.prediction
        # check that solution is indeed singular
        if !pull_back_is_to_affine(tracker.problem, state.solution)
            LinearAlgebra.normalize!(state.solution)
        end

        state.solution_cond = condition_jacobian(tracker)
        # If cond is not so high, maybe we don't have a singular solution in the end?
        if state.winding_number == 1 && state.solution_cond < 1e8
            result = correct!(core_tracker.state.x̄, core_tracker, state.solution, Inf;
                use_qr=true, max_iters=3,
                accuracy=core_tracker.options.refinement_accuracy)
            if isconverged(result)
                state.winding_number = 0
                @goto non_singular_case
            end
        end

        # In the case of a squared up system we now have to get rid of the
        # excess solutions, but since we don't have Newton's method at hand
        # we simply rely non the residual.
        if is_squared_up_system(core_tracker.homotopy) &&
           residual(tracker) > options.overdetermined_min_residual
             state.status = PathTrackerStatus.excess_solution
        end
    # 1.a) + 2.a)
    else
        # First we refine the obtained solution if possible
        result = correct!(core_tracker.state.x̄, core_tracker, currx(core_tracker), Inf;
            use_qr=true, max_iters=3,
            accuracy=core_tracker.options.refinement_accuracy)
        @label non_singular_case

        state.solution_cond = core_tracker.state.jacobian.cond

        if isconverged(result) && core_tracker.state.jacobian.corank_proposal == 0
            state.solution .= core_tracker.state.x̄
        elseif (!isconverged(result) || core_tracker.state.jacobian.corank_proposal > 0) &&
               at_infinity_post_check(state.val)

            state.solution .= currx(core_tracker)
            state.status = PathTrackerStatus.at_infinity
            return nothing
        else
            state.status = PathTrackerStatus.post_check_failed
            return nothing
        end

        # We want to check that if we present an affine solution to a user
        # that this is still in a convergent region of the affine target system
        # This covers that we have a squared up system and when we tracked in projective space
        if pull_back_is_to_affine(tracker.problem, state.solution) &&
            (is_squared_up_system(core_tracker.homotopy) ||
             !affine_tracking(core_tracker))

            # We start with bringing it on the affine patch
            if !affine_tracking(core_tracker)
                ProjectiveVectors.affine_chart!(state.solution)
            end
            init_auto_scaling!(cache.weighted_ip, state.solution, AutoScalingOptions())
            target_result = newton!(cache.base_point, cache.target_system, state.solution,
                            cache.weighted_ip, cache.target_newton_cache;
                            use_qr=true,
                            tol=core_tracker.options.refinement_accuracy, miniters=1, maxiters=2)
            if isconverged(target_result)
                state.solution .= cache.base_point
            elseif is_squared_up_system(core_tracker.homotopy)
                state.status = PathTrackerStatus.excess_solution
            else
                state.status = PathTrackerStatus.post_check_failed
            end
            state.solution_accuracy = target_result.accuracy

            if !pull_back_is_to_affine(tracker.problem, state.solution)
                LinearAlgebra.normalize!(state.solution)
                changepatch!(cache.target_system.patch, state.solution)
                state.solution_accuracy = result.accuracy
            else
                state.solution_accuracy = result.accuracy
            end
        # We have a purely projective solution. Here we only have to handle the case
        # of overdetermined systems.
        elseif !pull_back_is_to_affine(tracker.problem, state.solution) &&
                is_squared_up_system(core_tracker.homotopy)

            LinearAlgebra.normalize!(state.solution)
            changepatch!(cache.target_system.patch, state.solution)
            target_proj_result = newton!(cache.base_point, cache.target_system, state.solution,
                            euclidean_norm, cache.target_newton_cache;
                            use_qr=true,
                            tol=core_tracker.options.refinement_accuracy, miniters=1, maxiters=2)
            if isconverged(target_proj_result)
                state.solution .= cache.base_point
                state.solution_accuracy = result.accuracy
                LinearAlgebra.normalize!(state.solution)
            else
                state.status = PathTrackerStatus.excess_solution
            end

        end
    end
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
    tracker_startsolutions(prob, startsolutions; rest...)
end

function tracker_startsolutions(prob::Problem, startsolutions; kwargs...)
    core_tracker_supported, pathtracker_kwargs = splitkwargs(kwargs, coretracker_supported_keywords)
    tracker = PathTracker(prob, start_solution_sample(startsolutions); pathtracker_kwargs...)
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

    if return_code == :success
        t = 0.0
    else
        t = exp(-state.s)
    end
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
        valuation = copy(tracker.state.val.v)
        valuation_accuracy = abs.(tracker.state.val.v̇)
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
    jac = tracker.cache.target_jacobian
    jacobian!(jac.J, tracker.cache.target_system, x,
              tracker.cache.target_newton_cache.system_cache)
    updated_jacobian!(jac; update_infos=true)
    jac.cond
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
        if r.return_code != PathTrackerStatus.success
            println(io, " • t: $(r.t)")
        end
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
        if r.return_code != PathTrackerStatus.success
            println(io, " • t: $(r.t)")
        end
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
