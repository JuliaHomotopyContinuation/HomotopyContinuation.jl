# INTERFACE

# Each Strategy has it's own state and cache
abstract type MonodromyStrategy end
abstract type MonodromyStrategyParameters end
abstract type MonodromyStrategyCache end


"""
    parameters(strategy::MondromyStrategy, nparams::Integer)

Construct the parameters of the given `strategy`.
"""
function parameters end

"""
    regenerate(parameters::MonodromyStrategyParameters)::MonodromyStrategyParameters

Regenerate the parameters of given strategy.
"""
function regenerate end

"""
    cache(strategy::MonodromyStrategy, tracker)::MonodromyStrategyParameters

Regenerate the parameters of given strategy.
"""
function cache end

############
# Triangle
############
struct Triangle <: MonodromyStrategy end
struct TriangleParameters{N, T} <: MonodromyStrategyParameters
    p₁::SVector{N, T}
    p₂::SVector{N, T}

    γ::Union{Nothing, NTuple{2, ComplexF64}}
end

function parameters(strategy::Triangle, p₀::SVector{NParams, T}) where {NParams, T}
    p₁ = @SVector randn(T, NParams)
    p₂ = @SVector randn(T, NParams)
    γ = nothing
    if T <: Real
        γ = (randn(ComplexF64), randn(ComplexF64))
    end

    TriangleParameters(p₁, p₂, γ)
end

function regenerate(params::TriangleParameters{N, T}) where {N, T}
    p₁ = @SVector randn(T, N)
    p₂ = @SVector randn(T, N)
    γ = nothing
    if params.γ !== nothing
        γ = (randn(ComplexF64), randn(ComplexF64))
    end
    TriangleParameters(p₁, p₂, γ)
end

struct TriangleCache{T, H} <: MonodromyStrategyCache
    x₁::ProjectiveVectors.PVector{T, H}
    x₂::ProjectiveVectors.PVector{T, H}
end

function cache(strategy::Triangle, tracker::PathTracking.PathTracker)
    TriangleCache(copy(tracker.x), copy(tracker.x))
end


function loop(tracker, x₀::Vector, p₀::SVector, params::TriangleParameters, cache::TriangleCache, stats)
    H = Homotopies.basehomotopy(tracker.homotopy)::Homotopies.ParameterHomotopy

    Homotopies.set_parameters!(H, (p₀, params.p₁), params.γ)
    retcode = track!(cache.x₁, tracker, x₀, stats)
    if retcode != :success
        return x₀, retcode
    end

    Homotopies.set_parameters!(H, (params.p₁, params.p₂))
    retcode = track!(cache.x₁, tracker, cache.x₁, stats)
    if retcode != :success
        return x₀, retcode
    end

    Homotopies.set_parameters!(H, (params.p₂, p₀))
    retcode = track!(cache.x₁, tracker, cache.x₁, stats)

    return ProjectiveVectors.affine(cache.x₁), retcode
end

function track!(u, tracker::PathTracking.PathTracker, x, stats::MonodromyStatistics)
    retcode = PathTracking.track!(u, tracker, x, 1.0, 0.0)
    pathtracked!(stats, retcode)
    retcode
end
