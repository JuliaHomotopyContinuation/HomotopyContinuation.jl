# INTERFACE

# Each Strategy has it's own state and cache
abstract type AbstractStrategy end
abstract type AbstractStrategyParameters end
abstract type AbstractStrategyCache end

"""
    parameters(strategy::MondromyStrategy, nparams::Integer)

Construct the parameters of the given `strategy`.
"""
function parameters end

"""
    regenerate!(parameters::AbstractStrategyParameters)::AbstractStrategyParameters

Regenerate the parameters of given strategy.
"""
function regenerate end

"""
    cache(strategy::AbstractStrategy, tracker)::AbstractStrategyParameters

Regenerate the parameters of given strategy.
"""
function cache end


############
# Triangle
############
struct Triangle <: AbstractStrategy
    usegamma::Bool
end

Triangle(;usegamma=true) = Triangle(usegamma)

function graph(strategy::Triangle, p::SVector, x::AbstractVector, options::Options)
    Graph(p, x, 3, options, usegamma=strategy.usegamma)
end

struct DoubleEdge <: AbstractStrategy end

function graph(strategy::DoubleEdge, p::SVector, x::AbstractVector, options::Options)
    Graph(p, x, 2, options, usegamma=true)
end
