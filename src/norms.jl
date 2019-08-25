"""
    AbstractNorm

An `AbstractNorm` represents any norm of a vector space.
"""
abstract type AbstractNorm end

"""
    distance(u, v, norm::AbstractNorm)

Compute the distance ||u-v|| with the given norm.
"""
function distance end

"""
    WeightedNorm(norm, n::Integer)
    WeightedNorm(norm, x::AbstractVector)

A `WeightedNorm` represents a weighted variant of norm `norm` of a `n`-dimensional vector space.`
"""
struct WeightedNorm{N<:AbstractNorm} <: AbstractNorm
    weights::Vector{Float64}
    norm::N
end
WeightedNorm(norm::AbstractNorm, x::AbstractVector) = WeightedNorm(norm, length(x))
WeightedNorm(norm::AbstractNorm, n::Integer) = WeightedNorm(ones(n), norm)

(N::WeightedNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::WeightedNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    weights(WN::WeightedNorm)

Returns the inverse weights of the weighted norm.
"""
weights(WN::WeightedNorm) = WN.weights

Base.length(WN::WeightedNorm) = length(WN.weights)
Base.size(WN::WeightedNorm) = (length(WN),)
Base.copyto!(WN::WeightedNorm, x) = copyto!(WN.weights, x)
Base.@propagate_inbounds Base.getindex(WN::WeightedNorm, i::Integer) = getindex(WN.weights, i)
Base.@propagate_inbounds Base.setindex!(WN::WeightedNorm, x, i::Integer) = setindex!(WN.weights, x, i)

"""
    EuclideanNorm <: AbstractNorm

The usual euclidean norm resp. 2-norm.
"""
struct EuclideanNorm <: AbstractNorm end

function distance(x::AbstractVector, y::AbstractVector, ::EuclideanNorm)
    n = length(x)
    @boundscheck n == length(y)
    @inbounds d = abs2(x[1] - y[1])
    for i in 2:n
        @inbounds d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
function distance(x::AbstractVector, y::AbstractVector, w::WeightedNorm{EuclideanNorm})
    @boundscheck length(w) == length(x) == length(y)
    @inbounds d = abs2(x[1] - y[1]) / (w[1]^2)
    for i in 2:length(x)
        @inbounds d += abs2(x[i] - y[i]) / (w[i]^2)
    end
    sqrt(d)
end

function LinearAlgebra.norm(x::AbstractVector, ::EuclideanNorm)
    n = length(x)
    @inbounds d = abs2(x[1])
    for i in 2:n
        @inbounds d += abs2(x[i])
    end
    sqrt(d)
end
function LinearAlgebra.norm(x::AbstractVector, w::WeightedNorm{EuclideanNorm})
    @boundscheck length(w) == length(x)
    @inbounds out = abs2(x[1]) / (w[1]^2)
    for i in 2:length(x)
        @inbounds out += abs2(x[i]) / (w[i]^2)
    end
    sqrt(out)
end
(N::EuclideanNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::EuclideanNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    InfNorm <: AbstractNorm

The ∞-norm.
"""
struct InfNorm <: AbstractNorm end

function distance(x::AbstractVector, y::AbstractVector, ::InfNorm)
    n = length(x)
    @boundscheck n == length(y)
    @inbounds dmax = abs2(x[1] - y[1])
    for i in 2:n
        @inbounds dᵢ = abs2(x[i] - y[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
function distance(x::AbstractVector, y::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2(x[1] - y[1]) / (w[1]^2)
    for i in 2:n
        @inbounds dᵢ = abs2(x[i] - y[i]) / (w[i]^2)
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end

function LinearAlgebra.norm(x::AbstractVector, ::InfNorm)
    n = length(x)
    @inbounds dmax = abs2(x[1])
    for i in 2:n
        @inbounds dᵢ = abs2(x[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
function LinearAlgebra.norm(x::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2(x[1]) / (w[1]^2)
    for i in 2:n
        @inbounds dᵢ = abs2(x[i]) / (w[i]^2)
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
(N::InfNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::InfNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    euclidean_distance(u, v)

Compute ||u-v||₂.
"""
euclidean_distance(u, v) = distance(u, v, EuclideanNorm())

"""
    euclidean_norm(u)

Compute ||u||₂.
"""
euclidean_norm(x::AbstractVector) = LinearAlgebra.norm(x, EuclideanNorm())
