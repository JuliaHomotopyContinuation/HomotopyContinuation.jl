export NullPredictor

"""
    NullPredictor()

A predictor which does no prediction step, i.e., it just returns the input as
its prediction.
"""
struct NullPredictor <: AbstractPredictor end
struct NullPredictorCache <: AbstractPredictorCache end

cache(::NullPredictor, H, x, ẋ, t) = NullPredictorCache()

function predict!(xnext, ::NullPredictorCache, H, x, t, dt, ẋ, Jac)
    xnext .= x
    nothing
end

order(::NullPredictorCache) = 1
