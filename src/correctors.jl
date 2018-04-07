module Correctors

abstract type AbstractCorrector end
abstract type AbstractCorrectorCache end

"""
    cache(::AbstractCorrector, ::HomotopyWithCache{M, N}, x, t)::AbstractCorrectorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    correct!(xnext, ::AbstractCorrector, ::AbstractCorrectorCache, H::HomotopyWithCache, x, t)
"""
function correct! end

end
