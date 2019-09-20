"""
    CauchyEndgame(x::AbstractVector)

The Cauchy endgame tries to to predict the value a solution path `x(t)` at `0` , i.e. `x(0)`,
by using a generalization of [Cauchy's integral formula]. This is derived in [^MSW90].
The general idea is that we make loops around the origin until we have again the start point.
The number of loops is referred to as the *winding* or *cycle number*.

[^MSW90]: Morgan, A.P., Sommese, A.J. & Wampler, C.W. Numer. Math. (1990) 58: 669. https://doi.org/10.1007/BF01385648
[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
struct CauchyEndgame{AV<:AbstractVector}
    base_point::AV
    CauchyEndgame(x::AV) where {AV<:AbstractVector} = new{AV}(copy(x))
end

"""
    CauchyEndgameResult

An enum indicating the result of the [`predict!`](@ref) computation.

# Cases
* `CAUCHY_SUCCESS`: The endgame was successfull.
* `CAUCHY_TERMINATED_MAX_WINDING_NUMBER`: The endgame was terminated since the winding
  number is larger than the provided threshold.
* `CAUCHY_TERMINATED_ACCURACY_LIMIT`: The endgame was terminated since during the loop
  tracking the achievable accuracy was no more sufficient.
* `CAUCHY_TERMINATED`: The endgame was terminated due to some other error in the path
  tracking.
"""
@enum CauchyEndgameResult begin
    CAUCHY_SUCCESS
    CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    CAUCHY_TERMINATED_ACCURACY_LIMIT
    CAUCHY_TERMINATED
end

"""
    predict!(p, tracker::CoreTracker, ::CauchyEndgame;
             samples_per_loop::Int = 8, max_winding_number::Int = 12)

Try to predict the value of `x(0)` using the [`CauchyEndgame`](@ref).
For this we track the polygon defined by ``te^{i2πk/n}`` until we end again at ``x``.
Here ``n`` is the number of samples we take per loop, `samples_per_loop`.
The computation gives up if we have a winding number larger than `max_winding_number`.
It returns a tuple denoting the success ([`CauchyEndgameResult`](@ref)) and the computed
winding number `m`.

[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
function predict!(
    p,
    tracker::CoreTracker,
    eg::CauchyEndgame;
    samples_per_loop::Int = 8,
    max_winding_number::Int = 12,
)
    @unpack base_point = eg
    tracker.options.logarithmic_time_scale || throw(ArgumentError("The `CoreTracker` needs to use a logarithmic time scale for the cauchy endgame."))

    s = current_t(tracker)
    s_target = tracker.state.segment.target

    # store initial tracking information
    initial_segment = tracker.state.segment
    initial_s = tracker.state.s
    initial_Δs = tracker.state.Δs
    initial_step_size = tracker.options.initial_step_size


    # store the current value in p
    base_point .= current_x(tracker)
    p .= 0.0
    x = current_x(tracker)

    # during the loop we fix the affine patch
    fix_patch!(tracker)

    m = k = 1
    Δθ = 2π / samples_per_loop
    result = CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    while m ≤ max_winding_number
        θⱼ = 0.0
        for j = 1:samples_per_loop
            θⱼ₋₁ = θⱼ
            θⱼ += Δθ

            retcode = track!(tracker, x, s + im * θⱼ₋₁, s + im * θⱼ; loop = true)

            if !is_success(retcode)
                init!(tracker, base_point, s, s_target; keep_steps = true)
                if retcode == CT_TERMINATED_ACCURACY_LIMIT
                    return CAUCHY_TERMINATED_ACCURACY_LIMIT, m
                else
                    return CAUCHY_TERMINATED, m
                end
            end

            p .+= x
        end
        # Check that loop is closed
        if norm(tracker)(base_point, x) < 4 * tracker.options.accuracy
            p ./= m .* samples_per_loop
            result = CAUCHY_SUCCESS
            break
        end

        m += 1
    end

    # we have to undo the fixing of the patch
    unfix_patch!(tracker)
    init!(tracker, base_point, s, s_target; loop = true, keep_steps = true)

    result, m
end

fix_patch!(tracker::CoreTracker) = tracker.options.update_patch = false
unfix_patch!(tracker::CoreTracker) = tracker.options.update_patch = true
