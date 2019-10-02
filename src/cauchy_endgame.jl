"""
    CauchyEndgame(x::AbstractVector; samples_per_loop = 8, max_winding_number = 12)

The Cauchy endgame tries to to predict the value a solution path `x(t)` at `0` , i.e. `x(0)`,
by using a generalization of [Cauchy's integral formula]. This is derived in [^MSW90].
The general idea is that we make loops around the origin until we have again the start point.
The number of loops is referred to as the *winding* or *cycle number*.

[^MSW90]: Morgan, A.P., Sommese, A.J. & Wampler, C.W. Numer. Math. (1990) 58: 669. https://doi.org/10.1007/BF01385648
[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
struct CauchyEndgame{AV<:AbstractVector}
    base_point::AV
    samples_per_loop::Int
    max_winding_number::Int

    function CauchyEndgame(
        x::AV;
        samples_per_loop::Int = 8,
        max_winding_number::Int = 12,
    ) where {AV<:AbstractVector}
        new{AV}(copy(x), samples_per_loop, max_winding_number)
    end
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
    predict!(p, tracker::CoreTracker, ::CauchyEndgame)

Try to predict the value of `x(0)` using the [`CauchyEndgame`](@ref).
For this we track the polygon defined by ``te^{i2πk/n}`` until we end again at ``x``.
Here ``n`` is the number of samples we take per loop, `samples_per_loop`.
The computation gives up if we have a winding number larger than `max_winding_number`.
It returns a tuple denoting the success ([`CauchyEndgameResult`](@ref)) the computed
winding number `m::Int` and th expected accuracy of the solution.

[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
function predict!(p, tracker::CoreTracker, cauchy::CauchyEndgame)
    @unpack base_point, samples_per_loop, max_winding_number = cauchy

    tracker.options.logarithmic_time_scale || throw(ArgumentError("The `CoreTracker` needs to use a logarithmic time scale for the cauchy endgame."))

    s = current_t(tracker)
    s_target = tracker.state.segment.target

    # store the current value in p
    base_point .= current_x(tracker)
    p .= 0.0
    x = current_x(tracker)

    # during the loop we fix the affine patch
    fix_patch!(tracker)

    m = k = 1
    Δθ = 2π / samples_per_loop
    result = CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    # The accuracy of the Cauchy endgame depends on the minimal accuracy
    # of a sample point. If all sample points are computed with an accuracy τ
    # we can expect τ as an error of the solution
    max_acc = tracker.state.accuracy
    while m ≤ max_winding_number
        θⱼ = 0.0
        for j = 1:samples_per_loop
            θⱼ₋₁ = θⱼ
            θⱼ += Δθ

            retcode = track!(tracker, x, s + im * θⱼ₋₁, s + im * θⱼ; loop = true)
            max_acc = max(max_acc, tracker.state.accuracy)

            if !is_success(retcode)
                init!(tracker, base_point, s, s_target; keep_steps = true)
                if retcode == CoreTrackerStatus.terminated_accuracy_limit
                    return CAUCHY_TERMINATED_ACCURACY_LIMIT, m, max_acc
                else
                    return CAUCHY_TERMINATED, m, max_acc
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

    result, m, max_acc
end

fix_patch!(tracker::CoreTracker) = tracker.options.update_patch = false
unfix_patch!(tracker::CoreTracker) = tracker.options.update_patch = true
