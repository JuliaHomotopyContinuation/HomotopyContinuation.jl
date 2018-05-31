module PredictionCorrection

import ..Homotopies
import ..Predictors
import ..Correctors
using ..Utilities

export PredictorCorrector,
    PredictorCorrectorCache,
    cache,
    predict_correct!

"""
    PredictorCorrector(predictor::Predictors.AbstractPredictor, corrector::Correctors.AbstractCorrector)

Construct a predictor correction from the given `predictor` and `corrector`. This is mostly
a convenience wrapper to manage `predictor` and `corrector`.
"""
struct PredictorCorrector{P<:Predictors.AbstractPredictor, C<:Correctors.AbstractCorrector}
    predictor::P
    corrector::C
end

"""
    PredictorCorrector{P<:Predictors.AbstractPredictorCache, C<:Correctors.AbstractCorrectorCache}

A wrapper holding the predictor and corrector patch.
"""
struct PredictorCorrectorCache{V<:AbstractVector, P<:Predictors.AbstractPredictorCache, C<:Correctors.AbstractCorrectorCache}
    xnext::V
    predictor::P
    corrector::C
end

"""
    cache(PC::PredictorCorrector, H::HomotopyWithCache, x, t)

Assemble the cache for `PC`.
"""
function cache(PC::PredictorCorrector, H, x, t)
    PredictorCorrectorCache(
        copy(x),
        Predictors.cache(PC.predictor, H, x, t),
        Correctors.cache(PC.corrector, H, x, t)
    )
end

"""
    predict_correct!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H::HomotopyWithCache, t, Δt, tol, maxiters)

Perform a prediction-correction step and store the result in `x` if successfull. Returns `(true, :ok)` if
the correction was sucessfull, otherwise `(false, status)` where `status` is either `:ok` or an error code.
"""
function predict_correct!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H, t, Δt, tol, maxiters)
    try
        xnext = cache.xnext
        Predictors.predict!(xnext, PC.predictor, cache.predictor, H, x, t, Δt)
        result = Correctors.correct!(xnext, PC.corrector, cache.corrector, H, xnext, t + Δt, tol, maxiters)
        if result.converged
            x .= xnext
            return (true, :ok)
        else
            # we have to reset any patch updates
            Homotopies.update!(H, x, t)
            return (false, :ok)
        end
    catch err
        if err isa InterruptException
            throw(err)
        elseif !(err isa Base.LinAlg.SingularException)
            warn("Cached the following expression thrown from Predictor or Corrector:")
            warn(err)
        end
        # we have to reset any patch updates
        Homotopies.update!(H, x, t)
        return (false, :predictor_corrector_failed)
    end
end

"""
    refine!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H::HomotopyWithCache, t, tol, maxiters)

Perform a correction step and store the result in `x` if the new residual is better than the initial one.
Returns the achieved residual.
"""
function refine!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H, t, tol, maxiters, res₀)
    if res₀ < tol
        return res₀
    end
    try
        result = Correctors.correct!(cache.xnext, PC.corrector, cache.corrector, H, x, t, tol, maxiters)
        if result.res < res₀
            x .= cache.xnext
            return result.res
        end
        return res₀
    catch
        return res₀
    end
end

end
