module PatchSwitching

import ..AffinePatches
import ..Homotopies
import ..PathTracking
import ..Predictors
import ..StepLength
import ..Systems


"""
    PatchSwitcher(H, patch₁, patch₀, x, t; trackerkwargs...)

Construct a `PatchSwitcher` to switch the affine patch from `patch₁` to `patch₀`.
`H` is a homotopy such that ``H(x, t) ≈ 0`` and such that `x` is on the affine
patch defined by `patch₁`.
To actually switch patches use `[switch!](@ref)`.
"""
struct PatchSwitcher{T<:PathTracking.PathTracker, V<:AbstractVector}
    tracker::T
    start::V
end

function PatchSwitcher(H::Homotopies.AbstractHomotopy,
    p₁::AffinePatches.AbstractAffinePatchState,
    p₀::AffinePatches.AbstractAffinePatchState, x, t;
    steplength=StepLength.HeuristicStepLength(initial=1.0),
    predictor=Predictors.Euler(),
    kwargs...)
    fixed = Systems.FixedHomotopy(H, t)
    start = copy(x)
    homotopy = Homotopies.PatchSwitcherHomotopy(fixed, p₁, p₀)

    tracker = PathTracking.PathTracker(homotopy, x, 1.0, 0.0;
        steplength=steplength, predictor=predictor, kwargs...)
    PatchSwitcher(tracker, start)
end

"""
    switch!(x, patchswitcher)

Move `x` from it's current patch to the patch defined by `patchswitcher`.
Returns a `Symbol` indicating whether this was successfull. If so `:success`
is returned. Only if the switch was successfull `x` is updated with the new value.
"""
function switch!(x, switcher::PatchSwitcher)
    switcher.start .= x
    PathTracking.track!(x, switcher.tracker, switcher.start, 1.0, 0.0)
end


end
