export PathResult, PathTrackerStatus, PathTracker,
       pathtracker, pathtracker_startsolutions, solution,
       accuracy, residual, start_solution, is_success, is_failed, is_at_infinity,
       is_singular, is_nonsingular, is_real, is_projective, is_affine,
       set_parameters!, multiplicity


const pathtracker_supported_keywords = [
    :at_infinity_check, :samples_per_loop, :max_winding_number,
    :overdetermined_min_accuracy, :overdetermined_min_residual,
    :cond_eg_start, :min_cond_at_infinity,
    :t_eg_start, :tol_val_inf_accurate, :tol_val_finite_accurate, :accuracy, :accuracy_eg]

###############
## VALUATION ##
###############

mutable struct Valuation
    v::Vector{Float64}
    v̇::Vector{Float64}
    v̈::Vector{Float64}
    v_data::NTuple{3,Vector{Float64}}
    s_data::NTuple{3,Float64}
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
    v_data = (zeros(n), zeros(n), zeros(n))
    s_data = (NaN, NaN, NaN)
    Valuation(v, v̇, v̈, v_data, s_data)
end

function Base.show(io::IO, val::Valuation)
    println(io, typeof(val), ":")
    for name in [:v, :v̇, :v̈]
        println(io, " • ", name, " → ", getfield(val, name))
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

function reset!(val::Valuation)
    val.v .= NaN
    val.v̇ .= NaN
    val.v̈ .= NaN
    val.s_data = (NaN, NaN, NaN)
    val
end

function update!(val::Valuation, z::Vector, ż, s, Δs::Float64, ::Val)
    _update!(val, z, ż, s, Δs)
end
function update!(val::Valuation, z::PVector, ż, s, Δs::Float64, ::Val{false})
    _update!(val, z, ż, s, Δs)
end
function update!(val::Valuation, z::PVector, ż, s, Δs::Float64, ::Val{true})
    _update_affine!(val, z, ż, s, Δs)
end
function _update!(val::Valuation, z::AbstractVector, ż, s, Δs::Float64)
    @unpack v, v_data, s_data = val
    (v₂, v₁, v₀), (s₂, s₁, _) = v_data, s_data
    v .= v₀ .= ν.(z, ż)
    # cycle data around
    val.v_data = (v₀, v₂, v₁)
    val.s_data = (s, s₂, s₁)

    !isnan(s₁) && finite_differences!(val)

    val
end
function _update_affine!(val::Valuation, z::PVector, ż, s, Δs::Float64)
    @unpack v, v_data, s_data = val
    (v₂, v₁, v₀), (s₂, s₁, _) = v_data, s_data
    k = 1
    for (rⱼ, j) in ProjectiveVectors.dimension_indices_homvars(z)
        vⱼ = ν(z[j], ż[j])
        for i in rⱼ
            v[k] = v₀[k] =  ν(z[i], ż[i]) - vⱼ
            k += 1
        end
    end
    # cycle data around
    val.v_data = (v₀, v₂, v₁)
    val.s_data = (s, s₂, s₁)

    !isnan(s₁) && finite_differences!(val)

    val
end

function ν(z::Complex, ż::Complex)
    x, y = reim(z); ẋ, ẏ = reim(ż);
    -(x * ẋ + y * ẏ) / abs2(z)
end

"""
    finite_differences!(val::Valuation)

Use a finite difference scheme to approximate the first and second derivative
of the valuation map.

Since we have a non-uniform grid, we need a more elaborate difference scheme.
The implementation follows the formulas derived in [^BS05]

[^BS05]: Bowen, M. K., and Ronald Smith. "Derivative formulae and errors for non-uniformly spaced points." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 461.2059 (2005): 1975-1997.
"""
function finite_differences!(val::Valuation)
    @unpack v_data, s_data, v̇, v̈ = val
    v₃, v₂, v₁ = v_data
    s₃, s₂, s₁ = s_data
    s₃₁ = s₃ - s₁
    a₃, a₁ = s₃ - s₂, s₁ - s₂
    for i in eachindex(v̇)
        v̇[i] = (a₃*v₁[i])/(a₁*s₃₁) - ((a₁+a₃)*v₂[i])/(a₁*a₃) - (a₁*v₃[i])/(s₃₁*a₃)
        v̈[i] = -2v₁[i]/((a₁)*s₃₁) + 2v₂[i]/(a₁*a₃) + 2v₃[i]/(s₃₁*a₃)
    end
    nothing
end

@inline function update!(val::Valuation, core_tracker::CoreTracker; at_infinity_check::Bool=true)
    z = core_tracker.state.x
    ż = core_tracker.state.ẋ
    Δs = core_tracker.state.Δs_prev
    s = core_tracker.state.s
    if at_infinity_check
        update!(val, z, ż, s, Δs, Val(true))
    else
        update!(val, z, ż, s, Δs, Val(false))
    end
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
    samples_per_loop::Int
    max_winding_number::Int
    # The minimal residual a solution needs to have to be considered
    # a solution of the original system (only applied for singular solutions)
    overdetermined_min_residual::Float64
    overdetermined_min_accuracy::Float64
    # minimial condition number where the endgame starts
    cond_eg_start::Float64
    min_cond_at_infinity::Float64
    # maximal t where the endgame starts
    t_eg_start::Float64
    tol_val_inf_accurate::Float64
    tol_val_finite_accurate::Float64
    accuracy::Float64
    accuracy_eg::Float64
end

function PathTrackerOptions(prob::Problem;
            at_infinity_check=true,
            samples_per_loop::Int=12,
            max_winding_number::Int=12,
            overdetermined_min_residual::Float64=1e-3,
            overdetermined_min_accuracy::Float64=1e-4,
            cond_eg_start::Float64=1e4,
            min_cond_at_infinity::Float64=1e7,
            t_eg_start::Float64=0.1,
            tol_val_inf_accurate::Float64=1e-4,
            tol_val_finite_accurate::Float64=1e-3,
            accuracy::Float64=error("You have to set `accuracy`"),
            accuracy_eg::Float64=min(accuracy, 1e-5))

    PathTrackerOptions(at_infinity_check, samples_per_loop, max_winding_number,
                       overdetermined_min_residual,
                       overdetermined_min_accuracy, cond_eg_start,
                       min_cond_at_infinity,
                       t_eg_start,
                       tol_val_inf_accurate, tol_val_finite_accurate, accuracy, accuracy_eg)
end

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts

mutable struct PathTrackerState{V<:AbstractVector}
    status::PathTrackerStatus.states
    s::Float64
    prediction::V
    solution::V
    solution_accuracy::Float64
    solution_cond::Float64
    winding_number::Int
    val::Valuation
end

function PathTrackerState(x; at_infinity_check::Bool=true)
    status = PathTrackerStatus.tracking
    s = 0.0
    prediction = copy(x)
    solution = copy(x)
    solution_accuracy = solution_cond = NaN
    winding_number = 0
    val = Valuation(x; at_infinity_check=at_infinity_check)

    PathTrackerState(status, s, prediction, solution,
                    solution_accuracy, solution_cond,
                    winding_number, val)
end


function reset!(state::PathTrackerState)
    state.status = PathTrackerStatus.tracking
    state.s = 0.0
    state.prediction .= zero(eltype(state.prediction))
    state.solution .= zero(eltype(state.prediction))
    state.solution_accuracy = state.solution_cond = NaN
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

    x = current_x(core_tracker)
    if is_squared_up_system(core_tracker.homotopy)
        F = squared_up_system(core_tracker.homotopy).F
    else
        F = FixedHomotopy(core_tracker.homotopy, Inf)
    end
    if affine_tracking(core_tracker)
        target_system = F
    elseif pull_back_is_to_affine(prob)
        target_system = PatchedSystem(F, state(EmbeddingPatch(), x))
    else # result is projective
        target_system = PatchedSystem(F, state(OrthogonalPatch(), x))
    end
    target_newton_cache = newton_cache(target_system, x)

    res = evaluate(target_system, current_x(core_tracker), target_newton_cache.system_cache)
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
with an endgame, i.e., it can also deal with singular solutions as well as diverging paths.
We call a diverged path a path going to infinity.
By convention a path is always tracked from t₁ > 0 towards 0.
During the path tracking an approximation of the valuation of a Puiseux series expansion of the solution is computed.
This is used to decide whether a path is diverging.
To compute singular solutions Cauchy's integral formula is used. There you have to trace out loops around the solution.
The number of loops necessary to arrive back at the start point is called the *winding number*.
In order to construct a `PathTracker` it is recommended to use the [`pathtracker`](@ref) and
[`pathtracker_startsolutions`](@ref) helper functions.
With a `PathTracker` constructed you can track a single path using the [`track`](@ref) method.
The result of this will be a [`PathResult`](@ref).

## Keyword arguments
`PathTracker` is a wrapper around [`CoreTracker`](@ref)
and thus it is possible to set all options which are available for [`CoreTracker`](@ref).
There are the following `PathTracker` specific options:

### General endgame parameters
* `accuracy_eg::Float64=min(accuracy, 1e-5))`: It is possible to change the accuracy during the path tracking. Usually you want lower the accuracy.
* `cond_eg_start::Float64=1e4`: The endgame is only started if the condition of the Jacobian is larger than this threshold.
* `max_winding_number::Int=12`: This limits the maximal number of loops taken in applying Cauchy's formula.
* `min_cond_at_infinity::Float64=1e7`: A path is declared as going to infinity only if it's Jacobian is also larger than this threshold.
* `samples_per_loop::Int=12`: To compute singular solutions Cauchy's integral formula is used. The accuracy of the solutions increases with the number of samples per loop.
* `t_eg_start::Float64=0.1`: The endgame starts only if `t` is smaller than this threshold.
* `tol_val_inf_accurate::Float64=1e-4`: A valuation which would result in a path declared as going to infinity is only accepted if the estimated accuracy of the valuation is less than this threshold.
* `tol_val_finite_accurate::Float64=1e-3`: A valuation which would result in a proper solution is only accepted if the estimated accuracy of the valuation is less than this threshold. This is only affects solutions where the path has at some point near 0 a condition number larger than `cond_eg_start`.

### Overdetermined system specific
* `overdetermined_min_accuracy=1e-5`: The minimal accuracy a non-singular solution needs to have to be considered a solution of the original system.
* `overdetermined_min_residual=1e-3`: The minimal residual a singular solution needs to have to be considered a solution of the original system.
"""
struct PathTracker{V<:AbstractVector, Prob<:AbstractProblem, PTC<:PathTrackerCache, CT<:CoreTracker}
    problem::Prob
    core_tracker::CT
    state::PathTrackerState{V}
    options::PathTrackerOptions
    cache::PTC
end

function PathTracker(prob::AbstractProblem, x::AbstractVector{<:Number};
                at_infinity_check=default_at_infinity_check(prob),
                accuracy=default_accuracy(prob),
                min_step_size=1e-30, kwargs...)

    core_tracker_supported, optionskwargs = splitkwargs(kwargs, coretracker_supported_keywords)
    core_tracker = CoreTracker(prob, x;
                        log_transform=true, predictor=Pade21(),
                        min_step_size=min_step_size,
                        accuracy=accuracy,
                        core_tracker_supported...)
    state = PathTrackerState(core_tracker.state.x; at_infinity_check=at_infinity_check)
    options = PathTrackerOptions(prob; at_infinity_check=at_infinity_check,
                                  accuracy=accuracy, optionskwargs...)
    cache = PathTrackerCache(prob, core_tracker)
    PathTracker(prob, core_tracker, state, options, cache)
end

Base.show(io::IO, tracker::PathTracker) = print(io, "PathTracker")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathTracker) = x


seed(PT::PathTracker) = PT.problem.seed
default_at_infinity_check(prob::Problem{AffineTracking}) = true
default_at_infinity_check(prob::Problem{ProjectiveTracking}) = homvars(prob) !== nothing
default_accuracy(prob::Problem) = pull_back_is_to_affine(prob) ? 1e-6 : 1e-8

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
function track!(tracker::PathTracker, x₁, s₁=0.0, s₀=-log(tracker.core_tracker.options.min_step_size);
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
    set_accuracy!(core_tracker, options.accuracy)
    reset!(state)
    # Handle the case that the start value is already a solution
    if residual(tracker, core_tracker.state.x) < 1e-14
        state.status = PathTrackerStatus.success
        state.prediction .= core_tracker.state.x
        state.s = Inf
        # Set winding_number to 1 to get into the possible singular solution case in
        # `check_and_refine_solution!`
        state.winding_number = 1
    end

    s_eg_start = -log(options.t_eg_start)
    while state.status == PathTrackerStatus.tracking
        step!(core_tracker)
        state.s = real(current_t(core_tracker))
        check_terminated!(core_tracker)

        if core_tracker.state.status == CoreTrackerStatus.success
            state.status = PathTrackerStatus.success
            break
        elseif core_tracker.state.status != CoreTrackerStatus.tracking
            state.status = PathTrackerStatus.status(core_tracker.state.status)
            break
        end

        # We only care if we moved forward
        core_tracker.state.last_step_failed && continue
        # update valuation and associated data
        update!(state.val, core_tracker; at_infinity_check=tracker.options.at_infinity_check)
        verdict = judge(state.val, options)
        # If we are too early, we don't work on the endgame
        state.s > s_eg_start || verdict == VALUATION_FINITE || continue
        # If the condition number is too low we also do not care
        start_eg(core_tracker, options) || continue
        # possibly reduce desired accuracy
        if tracker.options.at_infinity_check &&
           verdict == VALUATION_AT_INFINIY &&
           start_infinity_cutoff(core_tracker, options)

            state.status = PathTrackerStatus.at_infinity
            break
        end

        if verdict == VALUATION_FINITE
            set_accuracy!(core_tracker, options.accuracy_eg)
            retcode = predict_with_cauchy_integral_method!(state, core_tracker, options, cache)
            if retcode == :success
                state.status = PathTrackerStatus.success
                break
            elseif retcode == :max_winding_number
                continue
            else # path tracker failed during the loops -> break
                state.status = PathTrackerStatus.tracker_failed
                break
            end
        end
    end

    check_and_refine_solution!(tracker)
    # We have to set the number of blas threads to the previous value
    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)
    state.status
end

function start_eg(core_tracker, options)
    C = max(cond(core_tracker), exp10(digits_lost(core_tracker)))
    C > options.cond_eg_start
end

function start_infinity_cutoff(core_tracker, options)
    C = max(cond(core_tracker), exp10(digits_lost(core_tracker)))
    C > options.min_cond_at_infinity
end

"""
    judge(val::Valuation, J::Jacobian)::ValuationVerdict

Judge the current valuation.
"""
function judge(val::Valuation, options::PathTrackerOptions)
    @unpack v, v̇, v̈ = val
    finite = true
    for i in eachindex(v)
        # A coordinate goes to infinity if its valuation is negative.
        # We argue that a valuation is "stable" if it's first and second
        # derivative are "small".
        if v[i] < -0.1 &&
            abs(v̇[i]) < options.tol_val_inf_accurate &&
            abs(v̈[i]) < options.tol_val_inf_accurate

            return VALUATION_AT_INFINIY
        end

        if v[i] < -options.tol_val_finite_accurate ||
            !(abs(v̇[i]) < options.tol_val_finite_accurate &&
              abs(v̈[i]) < options.tol_val_finite_accurate)

            finite = false
        end
    end
    finite && return VALUATION_FINITE

    VALUATION_INDECISIVE
end

function at_infinity_post_check(val::Valuation, options::PathTrackerOptions)
    for (vᵢ, v̇ᵢ) in zip(val.v, val.v̇)
        if vᵢ < -0.1 && abs(v̇ᵢ) < sqrt(options.tol_val_inf_accurate)
            return true
        end
    end
    false
end

"""
    track(tracker::PathTracker, x₁, t₁::Float64=1.0; path_number::Int=1, details::Symbol=:default, options...)::PathResult

Track the path with start solution `x₁` from `t₁` towards `t=0`. The `details` options controls
the level of details of the informations available in the [`PathResult`](@ref).

Possible values for the options are
* `accuracy::Float64`
* `max_corrector_iters::Int`
* `max_steps::Int`
* `start_parameters::AbstractVector`
* `target_parameters::AbstractVector`
"""
function track(tracker::PathTracker, x₁, t₁=0.0, t₀=-log(tracker.core_tracker.options.min_step_size); path_number::Int=1, details::Symbol=:default, kwargs...)
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

function set_eg_accuracy!(core_tracker::CoreTracker, options::PathTrackerOptions)
    set_accuracy!(core_tracker, options.accuracy_eg)
end


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

    initial_segment = core_tracker.state.segment
    initial_s = core_tracker.state.s
    initial_Δs = core_tracker.state.Δs
    initial_Δs_prev = core_tracker.state.Δs_prev
    initial_step_size = core_tracker.options.initial_step_size
    s = real(current_t(core_tracker))

    base_point .= current_x(core_tracker)
    prediction .= zero(eltype(state.prediction))

    # during the loop we fix the affine patch
    fix_patch!(core_tracker)

    m = k = 1
    ∂θ = 2π / samples_per_loop
    core_tracker.options.initial_step_size = ∂θ #0.5∂θ
    while m ≤ max_winding_number
        θⱼ = 0.0
        for j=1:samples_per_loop
            θⱼ₋₁ = θⱼ
            θⱼ += ∂θ

            retcode = track!(core_tracker, current_x(core_tracker), s + im*θⱼ₋₁, s + im*θⱼ; loop=true)
            if retcode != CoreTrackerStatus.success
                # during the loop we fixed the affine patch
                unfix_patch!(core_tracker)
                core_tracker.state.segment = initial_segment
                core_tracker.options.initial_step_size = initial_step_size
                core_tracker.state.s = initial_s
                core_tracker.state.Δs = initial_Δs
                core_tracker.state.Δs_prev = initial_Δs_prev
                core_tracker.state.status = CoreTrackerStatus.tracking
                return Symbol(retcode)
            end

            prediction .+= current_x(core_tracker)
        end

        if distance(base_point, current_x(core_tracker), inner(core_tracker)) < 4core_tracker.options.accuracy
            break
        end

        m += 1
    end
    core_tracker.state.segment = initial_segment
    core_tracker.options.initial_step_size = initial_step_size
    core_tracker.state.s = initial_s
    core_tracker.state.Δs = initial_Δs
    core_tracker.state.Δs_prev = initial_Δs_prev
    core_tracker.state.status = CoreTrackerStatus.tracking
    # we have to undo the fixing of the patch
    unfix_patch!(core_tracker)

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
        at_infinity_post_check(state.val, options)

       state.status = PathTrackerStatus.at_infinity
   end

    if state.status ≠ PathTrackerStatus.success
        state.solution .= current_x(core_tracker)
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
        if !pull_back_is_to_affine(tracker.problem)
            LinearAlgebra.normalize!(state.solution)
        end

        state.solution_cond = condition_jacobian(tracker)
        # If cond is not so high, maybe we don't have a singular solution in the end?
        if state.winding_number == 1 && state.solution_cond < 1e10
            try
                result = correct!(core_tracker.state.x̄, core_tracker, state.solution, Inf;
                    use_qr=true, max_iters=3,
                    precision=PRECISION_ADAPTIVE,
                    accuracy=core_tracker.options.refinement_accuracy)
            catch e
                result = correct!(core_tracker.state.x̄, core_tracker, state.solution, Inf;
                    use_qr=true, max_iters=3,
                    precision=PRECISION_FIXED_64,
                    accuracy=core_tracker.options.refinement_accuracy)
            end
            if isconverged(result)
                @goto non_singular_case
            end
        elseif state.winding_number > 1 && state.solution_cond < 1e10 && residual(tracker) > 100
            state.status = PathTrackerStatus.post_check_failed
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
        result = correct!(core_tracker.state.x̄, core_tracker, current_x(core_tracker), Inf;
            use_qr=true, max_iters=3,
            accuracy=core_tracker.options.refinement_accuracy)
        @label non_singular_case
        state.solution_cond = core_tracker.state.jacobian.cond

        if isconverged(result) && core_tracker.state.jacobian.corank_proposal == 0
            state.solution .= core_tracker.state.x̄
            state.solution_accuracy = result.accuracy
        elseif (!isconverged(result) || core_tracker.state.jacobian.corank_proposal > 0) &&
               at_infinity_post_check(state.val, options)

            state.solution .= current_x(core_tracker)
            state.status = PathTrackerStatus.at_infinity
            return nothing
        else
            state.solution .= current_x(core_tracker)
            state.status = PathTrackerStatus.post_check_failed
            return nothing
        end

        # We want to check that if we present an affine solution to a user
        # that this is still in a convergent region of the affine target system
        # This covers that we have a squared up system and when we tracked in projective space
        if pull_back_is_to_affine(tracker.problem) &&
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

            if !pull_back_is_to_affine(tracker.problem)
                LinearAlgebra.normalize!(state.solution)
                changepatch!(cache.target_system.patch, state.solution)
                state.solution_accuracy = result.accuracy
            else
                state.solution_accuracy = result.accuracy
            end
        # We have a purely projective solution. Here we only have to handle the case
        # of overdetermined systems.
        elseif !pull_back_is_to_affine(tracker.problem) &&
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
function pathtracker_startsolutions(args...; kwargs...)
    invalid = invalid_kwargs(kwargs, pathtracker_startsolutions_supported_keywords)
    check_kwargs_empty(invalid, pathtracker_startsolutions_supported_keywords)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    tracker_startsolutions(prob, startsolutions; rest...)
end

function tracker_startsolutions(prob::Problem, startsolutions; kwargs...)
    tracker = PathTracker(prob, start_solution_sample(startsolutions); kwargs...)
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
* `multiplicity::Union{Nothing, Int}` is the multiplicity of the `solution`. This is only assigned by. [`singular`](@ref).
* `condition_jacobian::Union{Nothing, Float64}`: This is the condition number of the row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.
* `winding_number:Union{Nothing, Int}`: The estimated winding number. This is a lower bound on the multiplicity of the solution.
* `start_solution::Union{Nothing, Int}`: The start solution of the path.
* `accepted_steps::Int`: The number of accepted steps during the path tracking.
* `rejected_steps::Int`: The number of rejected steps during the path tracking.
* `valuation::Union{Nothing, Vector{Float64}}`: An approximation of the valuation of the Puiseux series expansion of `x(t)`.
* `valuation_accuracy::Union{Nothing, Vector{Float64}}`: An estimate of the accuracy of the valuation of the Puiseux series expansion of `x(t)`.

     PathResult(tracker::PathTracker, start_solution=nothing, path_number::Union{Nothing,Int}=nothing; details=:default)

Possible `details` values are `:minimal` (minimal details), `:default` (default) and `:extensive` (all information possible).
"""
mutable struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Union{Nothing, Float64}
    residual::Union{Nothing, Float64} # level 1+
    multiplicity::Union{Nothing, Int} # only assigned by singular(Result)
    condition_jacobian::Union{Nothing, Float64}
    winding_number::Union{Nothing, Int}
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
        condition_jac = cond(core_tracker)
    else
        condition_jac = state.solution_cond
    end
    # residual
    if return_code == :success && details_level ≥ 1
        res = residual(tracker)
    else
        res = nothing
    end

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

    # this is only assigned by using the singular() function
    multiplicity = nothing

    PathResult(return_code, x, t, accuracy, res, multiplicity, condition_jac,
               windingnumber, path_number,
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
    multiplicity(P::PathResult{T})

Returns the multiplicity of `P`.
"""
multiplicity(P::PathResult{T}) where T = P.multiplicity

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
    is_success(pathresult)

Checks whether the path is successfull.
"""
is_success(r::PathResult) = r.return_code == :success

"""
    is_failed(pathresult)

Checks whether the path failed.
"""
is_failed(r::PathResult) =!(r.return_code == :at_infinity || r.return_code == :success)


"""
    is_at_infinity(pathresult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::PathResult) = r.return_code == :at_infinity

"""
    isfinite(pathresult)

Checks whether the path result is finite.
"""
Base.isfinite(r::PathResult) = r.return_code == :success # we don't check is_affine to make other code easier

"""
    is_singular(pathresult; tol=1e10)

Checks whether the path result is singular. This is true if
the multiplicity is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
is_singular(r::PathResult; tol=1e10) = is_singular(r, tol)
function is_singular(r::PathResult, tol::Real)
    (unpack(r.condition_jacobian, 1.0) > tol ||
     unpack(r.multiplicity, 1) > 1) &&
     is_success(r)
end

"""
    is_nonsingular(pathresult; tol=1e10)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
is_nonsingular(r::PathResult; kwargs...) = !is_singular(r; kwargs...) && is_success(r)
is_nonsingular(r::PathResult, tol::Real) = !is_singular(r, tol) && is_success(r)


"""
    is_real(pathresult; tol=1e-6)

We consider a result as `real` if the 2-norm of the imaginary part of the solution is at most `tol`.
"""
is_real(r::PathResult; tol=1e-6) = is_real(r, tol)
is_real(r::PathResult, tol::Real) = is_real_vector(r.solution, tol)
# provide fallback since this in in Base
Base.isreal(r::PathResult, tol) = is_real(r, tol)
Base.isreal(r::PathResult; kwargs...) = is_real(r; kwargs...)

"""
    is_projective(pathresult)

Return`s true if the solution is a projective vector.
"""
is_projective(r::PathResult{<:PVector}) = true
is_projective(r::PathResult) = false

"""
    is_affine(pathresult)

Return`s true if the solution is an affine vector.
"""
is_affine(r::PathResult) = !is_projective(r)
