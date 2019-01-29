export NullPredictor

"""
    NullPredictor()

A predictor which does no prediction step, i.e., it just returns the input as
its prediction.
"""
struct NullPredictor <: AbstractPredictor end
struct NullPredictorCache <: AbstractPredictorCache end

cache(::NullPredictor, H, x, ẋ, t) = NullPredictorCache()

function predict!(xnext, ::NullPredictor, ::NullPredictorCache, H, x, t, dt, ẋ)
    xnext .= x
    nothing
end

order(::NullPredictor) = 1
