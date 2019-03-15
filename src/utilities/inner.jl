export AbstractInnerProduct, EuclideanIP, WeightedIP, euclidean_distance, euclidean_norm,
       init_weight!, change_weight!, norm2, inner, distance


abstract type AbstractInnerProduct end

struct EuclideanIP <: AbstractInnerProduct end

"""
    WeightedIP(weight)

Represents a weighted euclidean inner product.
"""
struct WeightedIP{T<:Real} <: AbstractInnerProduct
    weight::Vector{T} # The applied weights are 1 / w[i]
end
WeightedIP(x::AbstractVector) = WeightedIP(zeros(real(eltype(x)), length(x)))

"""
    distance(u, v, inner_product=EuclideanIP)

Compute the distance ||u-v|| with the norm induced by `inner_product`.
If `inner_product === nothing` the euclidean scalar product is used.
"""
function distance(x::AbstractVector, y::AbstractVector, ::EuclideanIP=EuclideanIP())
    @boundscheck length(x) == length(y)
    n = length(x)
    @inbounds d = abs2(x[1] - y[1])
    @inbounds for i in 2:n
        @fastmath d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
function distance(x::AbstractVector, y::AbstractVector, in::WeightedIP)
    @boundscheck length(in.weight) == length(x) == length(y)
    w = in.weight
    @inbounds d = abs2(x[1] - y[1]) / (w[1]^2)
    for i in 2:length(x)
        @inbounds d += abs2(x[i] - y[i]) / (w[i]^2)
    end
    sqrt(d)
end


function norm2(x::AbstractVector, ::EuclideanIP=EuclideanIP())
    out = zero(real(eltype(x)))
    for i in eachindex(x)
        @inbounds out += abs2(x[i])
    end
    out
end
function norm2(x::AbstractVector, in::WeightedIP)
    @boundscheck length(in.weight) == length(x)
    out = zero(real(eltype(x)))
    w = in.weight
    for i in eachindex(x)
        @inbounds out += abs2(x[i]) / (w[i]^2)
    end
    out
end
LinearAlgebra.norm(x::AbstractVector, ip::AbstractInnerProduct) = sqrt(norm2(x, ip))

(ip::EuclideanIP)(x::AbstractVector) = LinearAlgebra.norm(x, ip)
(ip::WeightedIP)(x::AbstractVector) = LinearAlgebra.norm(x, ip)

"""
    euclidean_distance(u, v)

Compute ||u-v||₂.
"""
euclidean_distance(u, v) = distance(u, v, EuclideanIP())

"""
    euclidean_norm(u)

Compute ||u||₂.
"""
euclidean_norm(x::AbstractVector) = LinearAlgebra.norm(x, EuclideanIP())
