module Predictors

abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    cache(::AbstractPredictor, ::HomotopyWithCache{M, N}, x, t)::AbstractPredictorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    predict!(xnext, ::AbstractPredictor, ::AbstractPredictorCache, H::HomotopyWithCache, x, t)
"""
function predict! end



end
