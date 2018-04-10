module PredictionCorrection

import ..Predictors
import ..Correctors

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
    predict_correct!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H::HomotopyWithCache, t, dt, tol)

Perform a prediction-correction step and store the result in xnext. Returns `true` if
the correction was sucessfull, otherwise `false`.
"""
@inline function predict_correct!(x, PC::PredictorCorrector, cache::PredictorCorrectorCache, H, t, dt, tol)
    try
        xnext = cache.xnext
        Predictors.predict!(xnext, PC.predictor, cache.predictor, H, x, t, dt)
        successfull = Correctors.correct!(xnext, PC.corrector, cache.corrector, H, xnext, t - dt, tol)
        if successfull
            x .= xnext
            return (true, :ok)
        else
            return (false, :ok)
        end
    catch err
        warn(err)
        return (false, :predictor_corrector_failed)
    end
end

end
