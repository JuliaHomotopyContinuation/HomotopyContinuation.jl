import ..Predictors
import ..Correctors

export PredictorCorrector,
    PredictorCorrectorCache,
    cache,
    step!

"""
    PredictorCorrector(predictor::Predictors.AbstractPredictor, corrector::Correctors.AbstractCorrector)

Construct a predictor correction from the given `predictor` and `corrector`. This is mostly
a convenience wrapper to manage `predictor` and `corrector`.
"""
struct PredictorCorrector{P<:Predictors.AbstractPredictor, C<:Correctors.AbstractCorrector}
    predictor::P
    corrector::C
end

struct PredictorCorrectorCache{P<:Predictors.AbstractPredictorCache, C<:Correctors.AbstractCorrectorCache}
    predictor::P
    corrector::C
end

"""
    cache(PC::PredictorCorrector, H::HomotopyWithCache, x, t)

Assemble the cache for `PC`.
"""
function cache(PC::PredictorCorrector, H, x, t)
    PredictorCorrectorCache(
        Predictors.cache(PC.predictor, H, x, t),
        Correctors.cache(PC.corrector, H, x, t)
    )
end

"""
    step!(xnext, PC::PredictorCorrector, cache::PredictorCorrectorCache, H::HomotopyWithCache, x, t, dt, tol)

Perform a prediction-correction step and store the result in xnext. Returns `true` if
the correction was sucessfull, otherwise `false`.
"""
function step!(xnext, PC::PredictorCorrector, cache::PredictorCorrectorCache, H, x, t, dt, tol)
    Predictors.predict!(xnext, PC.predictor, cache.predictor, H, x, t, dt)
    Correctors.correct!(xnext, PC.corrector, cache.corrector, H, xnext, t, tol)
end
