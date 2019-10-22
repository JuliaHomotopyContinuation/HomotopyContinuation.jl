export MixPredictor

"""
    MixPredictor()

This uses a Padé-approximation of type (2,1) for prediction.
"""
struct MixPredictor{P1<:AbstractPredictor,P2<:AbstractPredictor} <: AbstractPredictor
    predictor1::P1
    predictor2::P2
end
struct MixPredictorCache{
    P1<:AbstractPredictorCache,
    P2<:AbstractPredictorCache,
} <: AbstractPredictorCache
    predictor1::P1
    predictor2::P2
    predictor1_active::Base.RefValue{Bool}
end

function cache(MP::MixPredictor, H, x, ẋ, t)
    p1 = cache(MP.predictor1, H, x, ẋ, t)
    p2 = cache(MP.predictor2, H, x, ẋ, t)
    MixPredictorCache(p1, p2, Ref(true))
end


function update!(
    cache::MixPredictorCache,
    H,
    x,
    ẋ,
    t,
    Jac::JacobianMonitor,
    eval_err::Float64,
)
    if cache.predictor1_active[]
        update!(cache.predictor1, H, x, ẋ, t, Jac, eval_err)
    else
        update!(cache.predictor2, H, x, ẋ, t, Jac, eval_err)
    end
    cache
end

function predict!(x̂, cache::MixPredictorCache, H::HomotopyWithCache, x, t, Δt, ẋ, J)
    if cache.predictor1_active[]
        predict!(x̂, cache.predictor1, H, x, t, Δt, ẋ, J)
    else
        predict!(x̂, cache.predictor2, H, x, t, Δt, ẋ, J)
    end
end


for f in [:order, :highest_derivative, :second_derivative, :third_derivative]
    @eval begin
        function $f(cache::MixPredictorCache)
            cache.predictor1_active[] ? $f(cache.predictor1) : $f(cache.predictor2)
        end
    end
end
