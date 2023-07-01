include("valuation.jl")


struct ProjectivePathValue
    z::Vector{ComplexF64}
    ż::Vector{ComplexF64}
    # Normalizator
    α::Vector{ComplexF64}
    # Scaling factor: Result of x * α' + 1
    λ::Base.RefValue{ComplexF64}
end

function ProjectivePathValue(tracker::Tracker)
    n = length(tracker.state.x)
    z = zeros(ComplexF64, n)
    ż = zeros(ComplexF64, n)
    α = randn(ComplexF64, n)

    ProjectivePathValue(z, ż, α, Ref(0.0 + 0.0im))
end

function update!(p::ProjectivePathValue, tracker::Tracker)
    λ = inv(1 + LA.dot(p.α, tracker.state.x))
    p.λ[] = λ

    p.z .= λ .* tracker.state.x

    assign_ẋ!(p.ż, tracker.predictor)
    λ² = λ^2
    p.ż .*= λ²

    p
end


Base.@kwdef struct PathSamplePoints
    samples::Vector{Tuple{Vector{ComplexF64},Float64}} = []
end

function init!(S::PathSamplePoints)
    empty!(S.samples)
end

function take_path_sample!(S::PathSamplePoints, tracker::Tracker)
    # TODO: We should refine only later this is wasteful assuming most of the time the samples are not needed
    refine_current_solution!(tracker)
    t = real(tracker.state.t)
    x = copy(tracker.state.x)
    push!(S.samples, (x, t))
end

Base.last(S::PathSamplePoints) = last(S.samples)
Base.length(S::PathSamplePoints) = length(S.samples)
Base.isempty(S::PathSamplePoints) = isempty(S.samples)

function extrapolate(S::PathSamplePoints; winding_number::Int, max_samples::Int = 10)
    n = length(S)
    x0, abs_err = Richardson.extrapolate(
        S.samples[max(1, n - max_samples + 1):n],
        power = inv(winding_number),
    )
    norm_x0 = LA.norm(x0)
    if norm_x0 > max(1e-8, abs_err)
        err = abs_err / norm_x0
    else
        err = abs_err
    end

    x0, err
end


module EndgameCode

@enum codes begin
    success
    tracking
    tracker_terminated
    path_diverging
end

end

is_tracking(code::EndgameCode.codes) = code == EndgameCode.tracking
is_success(code::EndgameCode.codes) = code == EndgameCode.success


Base.@kwdef mutable struct DivergingCheckState
    diverging_starts::Vector{Float64}
    diverging_abs_coords_at_start::Vector{Float64}
    diverging_conds_at_start::Vector{Float64}
end

function DivergingCheckState(n::Integer)
    DivergingCheckState(
        diverging_starts = zeros(Float64, n),
        diverging_abs_coords_at_start = zeros(Float64, n),
        diverging_conds_at_start = zeros(Float64, n),
    )
end

function init!(state::DivergingCheckState)
    state.diverging_starts .= NaN
    state.diverging_abs_coords_at_start .= NaN
    state.diverging_conds_at_start .= NaN
end


Base.@kwdef mutable struct EndgameState
    code::EndgameCode.codes
    valuation::Valuation
    in_endgame_zone::Bool = false
    endgame_zone_start_at::Float64 = NaN
    norm_at_endgame_zone_start::Float64 = NaN
    steps_eg::Int = 0
    target_system_normalized_norm_derivative::Float64
    singular_check_row_scaling::Vector{Float64}
    projective_path_value::ProjectivePathValue
    sample_points::PathSamplePoints = PathSamplePoints()
    diverging_check::DivergingCheckState

    # solution properties
    solution::Vector{ComplexF64}
    solution_accuracy::Float64 = NaN
    singular::Bool = false
    winding_number::Int = 0
    rcond::Float64 = NaN
end

function EndgameState(tracker::Tracker)
    (normalized_norm_derivative, coordinate_norms) =
        sample_homotopy_weyl_norm(tracker.homotopy, 0.0)

    EndgameState(
        code = EndgameCode.tracking,
        valuation = Valuation(tracker.state.x),
        target_system_normalized_norm_derivative = normalized_norm_derivative,
        singular_check_row_scaling = inv.(coordinate_norms),
        projective_path_value = ProjectivePathValue(tracker),
        solution = similar(tracker.state.x),
        diverging_check = DivergingCheckState(length(tracker.state.x)),
    )
end

Base.show(io::IO, S::EndgameState) = print_fieldnames(io, S)

function init!(state::EndgameState)
    state.code = EndgameCode.tracking
    state.in_endgame_zone = false
    state.endgame_zone_start_at = NaN
    state.norm_at_endgame_zone_start = NaN
    state.steps_eg = 0
    state.solution .= NaN
    state.solution_accuracy = NaN
    state.winding_number = 0
    state.rcond = NaN
    state.singular = false
    init!(state.valuation)
    init!(state.sample_points)
    init!(state.diverging_check)
end


function sample_homotopy_weyl_norm(H::AbstractHomotopy, t; iters = 25)
    m, n = size(H)
    u = TaylorVector{2}(zeros(ComplexF64, 2, m))
    u₀, u₁ = vectors(u)

    coordinate_norms = zeros(Float64, m)
    norm₀ = 0.0
    norm₁ = 0.0
    for _ = 1:iters
        x = LA.normalize!(randn(ComplexF64, n))
        taylor!(u, Val(1), H, x, (t, 1.0))
        norm₀ += LA.norm(u₀)
        norm₁ += LA.norm(u₁)
        coordinate_norms .+= abs.(u₀)
    end
    norm₀ /= iters
    norm₁ /= iters
    normalized_norm_derivative = norm₁ / norm₀

    coordinate_norms ./= iters

    normalized_norm_derivative, coordinate_norms
end

"""
    AbstractPathTracker

Supertype for path trackers.
"""
abstract type AbstractPathTracker end

Base.broadcastable(T::AbstractPathTracker) = Ref(T)

"""
    EndgameOptions(; options...)

Options controlling the behaviour of a [`EndgameTracker`](@ref).

## Options

* `at_infinity_check = true`: Whether divering paths should be truncated.
* `endgame_start = 0.1`: The point `t` in time where the endgame starts. Set it to `0.0`
  to disable the endgame.
* `only_nonsingular = false`: If `true` don't run the Cauchy endgame to handle singular
  solutions.
* `zero_is_at_infinity = false`: Whether paths going to a solution where at least one
  coordinates is zero should also be considered diverging.

### Parameters
These parameters control the behaviour during the endgame.

* `max_winding_number = 6`: The maximal winding number which is considered to compute a singular soution.
* `singular_rcond_threshold = 1e-12`: Value for the inverse of the condition number under which a solution is considered singular
* `valuation_accuracy_tol = 0.01`: Tolerance on the accuracy of the valuation
"""
Base.@kwdef mutable struct EndgameOptions
    # options
    zero_is_at_infinity::Bool = false
    at_infinity_check::Bool = true
    only_nonsingular::Bool = false

    # parameters
    singular_rcond_threshold::Float64 = 1e-12
    diverging_rcond_threshold::Float64 = 1e-8
    min_norm_diverging_growth::Float64 = 1e4
    valuation_accuracy_tol::Float64 = 0.01
    max_winding_number::Int = 6
end


struct EndgameTracker{H<:AbstractHomotopy}
    tracker::Tracker{H}
    state::EndgameState
    options::EndgameOptions
end

EndgameTracker(tracker::Tracker, options::EndgameOptions) =
    EndgameTracker(tracker, EndgameState(tracker), options)
EndgameTracker(tracker::Tracker; options::EndgameOptions = EndgameOptions()) =
    EndgameTracker(tracker, options)

EndgameTracker(
    H::AbstractHomotopy;
    tracker_options::TrackerOptions = TrackerOptions(),
    options::EndgameOptions = EndgameOptions(),
) = EndgameTracker(Tracker(H, options = tracker_options), options)


Base.broadcastable(T::EndgameTracker) = Ref(T)

function check_in_endgame_zone(eg::EndgameTracker; tol::Float64 = 1e-2)
    tracker = eg.tracker
    t = real(tracker.state.t)

    # 1) Distance in the solution variety
    # ||(x(t), H(⋅,t)) - (x(0), H(⋅,0))|| ≈ t * (||x'(t)|| + ||H_t(⋅,0)||)
    dist =
        t * (
            norm(eg.state.projective_path_value.ż) +
            eg.state.target_system_normalized_norm_derivative
        )

    if dist < tol
        eg.state.in_endgame_zone = true
        if eg.state.steps_eg == 0
            eg.state.endgame_zone_start_at = t
            eg.state.norm_at_endgame_zone_start = norm(tracker.state.x, Inf)
        end
        eg.state.steps_eg += 1
    else
        eg.state.in_endgame_zone = false
    end

    eg.state.in_endgame_zone
end

function check_diverging(eg::EndgameTracker, val_result::ValuationAnalyzeResult)
    val_result.verdict == ValuationVerdict.Diverging || return false

    for i in eachindex(val_result.coordinate_verdicts)
        if check_diverging(eg, val_result, i)
            eg.state.code = EndgameCode.path_diverging
            return true
        end
    end

    false
end
function check_diverging(eg::EndgameTracker, val_result::ValuationAnalyzeResult, i::Integer)
    @unpack tracker, state, options = eg

    verdict = val_result.coordinate_verdicts[i]

    if verdict != ValuationVerdict.Diverging
        state.diverging_check.diverging_starts[i] = NaN
        return false
    end

    if isnan(state.diverging_check.diverging_starts[i])
        t = real(tracker.state.t)
        state.diverging_check.diverging_starts[i] = t
        state.diverging_check.diverging_abs_coords_at_start[i] =
            fast_abs(tracker.state.x[i])
        return false
    end

    xᵢ = tracker.state.x[i]

    coordinate_growth =
        fast_abs(xᵢ) / state.diverging_check.diverging_abs_coords_at_start[i]

    # Check if coordinate growth is larger than threshold
    valᵢ = state.valuation.val_x[i]
    # If a path is diverging and we are truly in the zone where the expansion as a Puiseux series holds
    # then it holds that x(t) / x(t_start) ≈ (t/t_start)^valᵢ [dividing removes the coefficient] 
    # If we require that coordinate growth holds for a period τ where τ = t/t_start, then we have
    # that coordinate_growth ≈ τ^valᵢ

    growth_base = options.min_norm_diverging_growth
    min_growth = 32
    min_growth_extended_prec = 2
    base_growth_threshold = inv(growth_base)^valᵢ
    growth_threshold = max(min_growth, base_growth_threshold)
    growth_threshold_extended_prec = max(min_growth_extended_prec, base_growth_threshold)

    if coordinate_growth > growth_threshold ||
       (tracker.state.extended_prec && coordinate_growth > growth_threshold_extended_prec)
        return true
    end

    false
end

function check_singular(eg::EndgameTracker, val_result::ValuationAnalyzeResult)
    @unpack tracker, state, options = eg

    (is_finite(val_result) && eg.tracker.state.extended_prec) || return false

    m, _ = estimate_winding_number(eg.state.valuation; max_winding_number = 6)

    if m > 1 &&
       check_valuation_accuracy(eg.state.valuation, m) &&
       length(eg.state.sample_points) >= 2m

        _, err_1 = extrapolate(eg.state.sample_points; winding_number = 1, max_samples = 3m)
        x0_m, err_m =
            extrapolate(eg.state.sample_points; winding_number = m, max_samples = 3m)
        # If m is correct than err_m should be meaningful smaller than err_1
        # We have the rough estimate err_m = t^{N/m} whereas otherwise err_1 = t^{1/m}
        if err_m < max(err_1^2, err_1 / 100)
            eg.state.solution .= x0_m
            eg.state.solution_accuracy = err_m
            eg.state.code = EndgameCode.success
            eg.state.singular = true
            eg.state.winding_number = m
            eg.state.rcond = compute_endpoint_rcond(eg, x0_m)
            return
        end

        # See if we can find a singular solution with winding number 1 - this is
        # the case if solution is on a positive dimensional component
        # In this case val(ẋ[i]) ≈ 0 for all i
    elseif m == 1 &&
           (eg.tracker.state.statistics.failed_attempts_to_go_to_zero >= 3) &&
           length(eg.state.sample_points) >= 3 &&
           check_valuation_regular_finite(eg.state.valuation)

        x0_1, err_1 = extrapolate(eg.state.sample_points; winding_number = 1)
        rcond = compute_endpoint_rcond(eg, x0_1)
        # Check if we have really a singular solution
        if rcond < eg.options.singular_rcond_threshold
            eg.state.solution .= s
            eg.state.solution_accuracy = err_1
            eg.state.code = EndgameCode.success_singular
            eg.state.singular = true
            eg.state.winding_number = 1
            eg.state.rcond = rcond
        end
    end

end

function endgame(eg::EndgameTracker)
    @unpack tracker, state, options = eg

    update!(state.projective_path_value, tracker)
    check_in_endgame_zone(eg) || return

    # Check if we moved forward
    !tracker.state.last_step_failed || return

    # Take sample point and update valuation
    take_path_sample!(eg.state.sample_points, tracker)
    update!(state.valuation, tracker.predictor, real(tracker.state.t))

    val_result = analyze_valuation(eg.state.valuation; zero_is_finite = true)

    check_diverging(eg, val_result)

    check_singular(eg, val_result)
end

function compute_endpoint_rcond(eg::EndgameTracker, x = eg.tracker.state.x)
    ws = eg.tracker.ws
    J = ws.A
    evaluate_and_jacobian!(ws.work_n, J, eg.tracker.homotopy, x, 0.0)
    # apply row and column scaling to J
    C = eg.tracker.state.norm.weights
    R = eg.state.singular_check_row_scaling
    m, n = size(J)
    @inbounds for j = 1:n
        cj = C[j]
        for i = 1:m
            J[i, j] = cj * R[i] * J[i, j]
        end
    end

    norm_J = norm(J, Inf)
    if norm_J < n * 1e-12
        rcond = norm_J
    else
        # disable automatic scaling temporarily since we control row scaling
        equi_option = ws.equilibriate[]
        ws.equilibriate[] = false
        set_A!(ws, J)
        rcond = rcond!(eg.tracker.ws)
        ws.equilibriate[] = equi_option
    end

    rcond
end

function handle_tracker_tracking_ended(eg::EndgameTracker)
    @unpack tracker = eg
    if is_success(status(tracker))
        eg.state.solution .= tracker.state.x
        eg.state.solution_accuracy = tracker.state.accuracy
        eg.state.winding_number = 1
        eg.state.code = EndgameCode.success
        rcond = compute_endpoint_rcond(eg)
        eg.state.rcond = rcond
        eg.state.singular = rcond < eg.options.singular_rcond_threshold
    else
        eg.state.code = EndgameCode.tracker_terminated
    end
end


function init!(
    eg::EndgameTracker,
    x::AbstractVector,
    t₁ = 1.0,
    t₀ = 0.0;
    keep_steps::Bool = false,
)
    init!(eg.tracker, x, t₁, t₀; keep_steps = keep_steps)
    init!(eg.state)
    eg
end

function step!(eg::EndgameTracker)
    step_success = step!(eg.tracker)

    if is_tracking(eg.tracker.state.code)
        endgame(eg)
    else
        handle_tracker_tracking_ended(eg)
    end

    step_success
end

function track!(
    eg::EndgameTracker,
    x::AbstractVector,
    t₁ = 1.0,
    t₀ = 0.0;
    keep_steps::Bool = false,
)
    init!(eg, x, t₁, t₀; keep_steps = keep_steps)

    while is_tracking(eg.state.code)
        step!(eg)
    end

    eg.state.code
end

function track(
    eg::EndgameTracker,
    x::AbstractVector,
    t₁ = 1.0;
    path_number::Union{Int,Nothing} = nothing,
)
    track!(eg, x, t₁)
    PathResult(eg, x, path_number)
end

function PathResult(
    endgame_tracker::EndgameTracker,
    start_solution = nothing,
    path_number = nothing,
)
    @unpack tracker, state, options = endgame_tracker
    H = tracker.homotopy

    if is_success(state.code)
        solution = get_solution(H, state.solution, 0.0)
        t = 0.0
        evaluate!(tracker.corrector.r, H, state.solution, complex(0.0))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    else
        t = real(tracker.state.t)
        solution = get_solution(H, tracker.state.x, t)
        evaluate!(tracker.corrector.r, H, tracker.state.x, complex(0.0))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    end

    last_path_point = nothing
    if !isempty(state.sample_points)
        last_path_point = last(state.sample_points)
    end


    PathResult(
        return_code = Symbol(state.code),
        solution = solution,
        t = t,
        singular = state.singular,
        accuracy = state.solution_accuracy,
        residual = residual,
        condition_jacobian = 1 / state.rcond,
        winding_number = state.winding_number,
        last_path_point = last_path_point,
        valuation = state.in_endgame_zone ? copy(state.valuation.val_x) : nothing,
        start_solution = start_solution,
        path_number = path_number,
        ω = tracker.state.ω,
        μ = tracker.state.μ,
        extended_precision = tracker.state.extended_prec,
        accepted_steps = tracker.state.statistics.accepted_steps,
        rejected_steps = tracker.state.statistics.rejected_steps,
        steps_eg = state.steps_eg,
        extended_precision_used = tracker.state.used_extended_prec,
    )
end


# EndgamePathIterator #
struct EndgamePathIterator{T<:EndgameTracker}
    tracker::T
    t_real::Bool
end
Base.IteratorSize(::Type{<:EndgamePathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:EndgamePathIterator}) = Base.HasEltype()

"""
    iterator(tracker::Tracker, x₁, t₁=1.0, t₀=0.0)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific
path. In each iteration the tuple `(x,t)` is returned.

## Example

Assume you have `Tracker` `tracker` and you wan to track `x₁` from 1.0 to 0.25:
```julia
for (x,t) in iterator(tracker, x₁, 1.0, 0.25)
    println("x at t=\$t:")
    println(x)
end
```

Note that this is a stateful iterator. You can still introspect the state of the tracker.
For example to check whether the tracker was successfull
(and did not terminate early due to some problem) you can do
```julia
println("Success: ", is_success(status(tracker)))
```
"""
function iterator(tracker::EndgameTracker, x₁, t₁ = 1.0, t₀ = 0.0; kwargs...)
    init!(tracker, x₁, t₁, t₀; kwargs...)
    EndgamePathIterator(tracker, typeof(t₁ - t₀) <: Real)
end

function current_x_t(iter::EndgamePathIterator)
    @unpack x, t = iter.tracker.tracker.state
    (copy(x), iter.t_real ? real(t) : t)
end

function Base.iterate(iter::EndgamePathIterator, state = nothing)
    state === nothing && return current_x_t(iter), 1
    !is_tracking(iter.tracker.state.code) && return nothing

    while is_tracking(iter.tracker.state.code)
        step_failed = !step!(iter.tracker)
        step_failed || break
    end
    current_x_t(iter), state + 1
end

include("path_info.jl")
