export InnerProduct, EuclideanInnerProduct,
       init_weight!, change_weight!, norm2, inner, distance

import ProjectiveVectors
import ProjectiveVectors: PVector
#################
# Inner product
#################

abstract type InnerProduct end

"""
    EuclideanInnerProduct(weight)

Represents a weighted euclidean inner product.
"""
struct EuclideanInnerProduct{T<:Real} <: InnerProduct
    weight::Vector{T} # The applied weights are 1 / w[i]
end

Base.@propagate_inbounds weight(IP::EuclideanInnerProduct, i::Integer) = IP.weight[i]
Base.@propagate_inbounds Base.getindex(IP::EuclideanInnerProduct, i) = weight(IP, i)

const SCALE_MIN = 0.01
const SCALE_ABS_MIN = min(SCALE_MIN^2, 200 * sqrt(eps()))
const SCALE_MAX = 1.0 / eps() / sqrt(2)

function init_weight!(ip::EuclideanInnerProduct, x::PVector{T,1}) where {T}
    point_norm = norm(x, ip)[1]
    w = ip.weight
    for i in 1:length(w)
        wᵢ = abs(x[i])
        if wᵢ < SCALE_MIN * point_norm
            wᵢ = max(SCALE_MIN * point_norm, SCALE_ABS_MIN)
        elseif wᵢ > SCALE_MAX * point_norm
            wᵢ = SCALE_MAX * point_norm
        end
        w[i] = wᵢ
    end
end

function change_weight!(ip::EuclideanInnerProduct, x::PVector{T,1}) where {T}
    point_norm = norm(x, ip)[1]
    @show point_norm
    w = ip.weight
    for i in 1:length(w)
        wᵢ = (abs(x[i]) + w[i]) / 2
        if wᵢ < SCALE_MIN * point_norm
            wᵢ = max(SCALE_MIN * point_norm, SCALE_ABS_MIN)
        elseif wᵢ > SCALE_MAX * point_norm
            wᵢ = SCALE_MAX * point_norm
        end
        w[i] = wᵢ
    end
end


function inner(in::EuclideanInnerProduct, x::AbstractVector, y::AbstractVector)
    @boundscheck length(in.weight) == length(x) == length(y)
    @inbounds inner(in, x, y, Base.OneTo(length(x)))
end

Base.@propagate_inbounds function inner(ip::EuclideanInnerProduct, x::AbstractVector, y::AbstractVector, range)
    n = length(x)
    w = ip.weight
    out = zero(promote_type(eltype(x), eltype(y)))
    for i in range
        out += x[i] * conj(y[i]) / (w[i]^2)
    end
    out
end
(ip::EuclideanInnerProduct)(x::AbstractVector, y::AbstractVector) = inner(ip, x, y)


function distance(x::AbstractVector, y::AbstractVector, in::EuclideanInnerProduct)
    @boundscheck length(in.weight) == length(x) == length(y)
    out = zero(real(eltype(x)))
    w = in.weight
    for i in eachindex(x)
        @inbounds out += abs2(x[i] - y[i]) / (w[i]^2)
    end
    sqrt(out)
end



function norm2(x::AbstractVector, in::EuclideanInnerProduct)
    @boundscheck length(in.weight) == length(x)
    out = zero(real(eltype(x)))
    w = in.weight
    for i in eachindex(x)
        @inbounds out += abs2(x[i]) / (w[i]^2)
    end
    out
end
function norm2(z::PVector{T, 1}, ip::EuclideanInnerProduct) where {T}
    (norm2(z.data, ip),)
end
@generated function norm2(z::PVector{T, N}, ip::EuclideanInnerProduct) where {T, N}
    quote
        r = ProjectiveVectors.dimension_indices(z)
        @inbounds $(Expr(:tuple, (:(_norm2_range(z, r[$i], ip)) for i=1:N)...))
    end
end
Base.@propagate_inbounds function _norm2_range(z::PVector{T}, rᵢ::UnitRange{Int}, ip::EuclideanInnerProduct) where {T}
    w = ip.weight
    normᵢ = zero(real(T))
    for k in rᵢ
        normᵢ += abs2(z[k]) / (w[k]^2)
    end
    normᵢ
end

LinearAlgebra.norm(x::AbstractVector, ip::EuclideanInnerProduct) = sqrt(norm2(x, ip))
LinearAlgebra.norm(x::PVector, ip::EuclideanInnerProduct) = sqrt.(norm2(x, ip))
