module ProjectiveVectors

using LinearAlgebra
import Base: ==
import ..Utilities: infinity_norm, unsafe_infinity_norm

export AbstractProjectiveVector,
    PVector,
    ProdPVector,
    raw,
    homvar,
    affine,
    affine!,
    embed,
    at_infinity,
    pvectors,
    infinity_norm,
    infinity_norm_fast,
    unsafe_infinity_norm,
    converteltype

abstract type AbstractProjectiveVector{T, N} <: AbstractVector{T} end

"""
    PVector(z, homvar)

A `PVector` represents a projective vector `z` which is the result
of the embedding of an affine vector `x` into projective space.
Note that the constructor *does not* perform the embedding but rather
is a simple wrapper. In order to embed an affine vector use
[`embed`](@ref).
"""
struct PVector{T, H<:Union{Nothing, Int}} <: AbstractProjectiveVector{T, 1}
    data::Vector{T}
    homvar::H

    function PVector{T, Int}(data, homvar::Int) where {T}
        @assert 1 ≤ homvar ≤ length(data)
        new(data, homvar)
    end

    function PVector{T, Nothing}(data, homvar::Nothing) where {T}
        new(data, nothing)
    end
end
function PVector(x::Vector{T}, homvar::H=nothing) where {T, H<:Union{Nothing, Int}}
    PVector{T, H}(x, homvar)
end

Base.copy(z::PVector) = PVector(copy(z.data), z.homvar)

"""
    raw(z::PVector)

access_input the underlying vector of `z`. This is useful to pass
the vector into some function which does not know the
projective structure.
"""
raw(z::PVector) = z.data
raw(z) = z

"""
    homvar(z::PVector)

Get the index of the homogenous variable.
"""
homvar(z::PVector) = z.homvar

(==)(v::PVector, w::PVector) = v.homvar == w.homvar && v.data == w.data

"""
    converteltype(v::PVector, T)
"""
Base.similar(v::PVector, ::Type{T}) where T = PVector(convert.(T, v.data), v.homvar)

"""
    embed(x::Vector, homvar::Union{Nothing, Int})::PVector

Embed a vector `x` into projective space with homogenization variable `homvar`.

    embed(z::PVector, x::Vector)

Embed a vector `x` into projective space the same way `x` was embedded.
"""
function embed(x::Vector{T}, homvar) where T
    k = 1
    data = Vector{T}(undef, length(x) + 1)
    @inbounds for k in 1:length(data)
        if k == homvar
            data[k] = one(T)
        else
            i = k < homvar ? k : k - 1
            data[k] = x[i]
        end
    end
    PVector(data, homvar)
end
embed(z::PVector, x::Vector) = embed(x, z.homvar)
# TODO define inplace variant of embed


"""
    affine_normalizer(z::PVector)

Returns a scalar to bring `z` on its affine patch.
"""
affine_normalizer(z::PVector{T, Int}) where T = inv(z.data[z.homvar])
"""
    abs2_ffine_normalizer(z::PVector)

Returns the squared absolute value of the normalizer.
"""
abs2_affine_normalizer(z::PVector{T, Int}) where T = inv(abs2(z.data[z.homvar]))

"""
    affine(z::PVector{T, Int})::Vector

Return the corresponding affine vector with respect to the standard patch.

    affine(z::PVector, i::Int)::Vector

Return the corresponding affine vector with xᵢ=1.
"""
affine(z::PVector{T, Int}) where {T} = affine(z, z.homvar)
function affine(z::PVector{T}, i::Int) where T
    x = Vector{T}(undef, length(z) - 1)
    normalizer = inv(z.data[i])
    @inbounds for k in 1:length(z)
        if k == i
            continue
        end
        j = k < i ? k : k - 1
        x[j] = z.data[k] * normalizer
    end
    x
end

"""
    affine!(z::PVector{T, Int})

Bring the projective vector on associated affine patch.
"""
affine!(z::PVector{T, Int}) where T = LinearAlgebra.rmul!(z.data, affine_normalizer(z))

"""
    infinity_norm(z::PVector{T, Int})

Compute the ∞-norm of `z`. If `z` is a complex vector this is more efficient
than `norm(z, Inf)`.

    infinity_norm(z₁::PVector{T, Int}, z₂::PVector{T, Int})

Compute the ∞-norm of `z₁-z₂` by bringing both vectors first on their respective
affine patch. This therefore only makes sense if both vectors have the
same affine patch.
"""
function infinity_norm(z::PVector{<:Complex, Int})
    sqrt(maximum(abs2, raw(z)) * abs2_affine_normalizer(z))
end

function infinity_norm(z₁::PVector{<:Complex, Int}, z₂::PVector{<:Complex, Int})
    normalizer₁ = affine_normalizer(z₁)
    normalizer₂ = -affine_normalizer(z₂)
    @inbounds m = abs2(muladd(z₁[1], normalizer₁, z₂[1] * normalizer₂))
    n₁, n₂ = length(z₁), length(z₂)
    if n₁ ≠ n₂
        return convert(typeof(m), Inf)
    end
    @inbounds for k=2:n₁
        m = max(m, abs2(muladd(z₁[k], normalizer₁, z₂[k] * normalizer₂)))
    end
    sqrt(m)
end

"""
    unsafe_infinity_norm(z₁::PVector, z₂::PVector)

Compute the ∞-norm of `z₁-z₂`. This *does not* bring both vectors on their respective
affine patch before computing the norm.
"""
function unsafe_infinity_norm(z₁::PVector, z₂::PVector)
    @inbounds m = abs2(z₁[1] - z₂[1])
    n₁, n₂ = length(z₁), length(z₂)
    if n₁ ≠ n₂
        return convert(typeof(m), Inf)
    end
    @inbounds for k=2:n₁
        m = max(m, abs2(z₁[k] - z₂[k]))
    end
    sqrt(m)
end

"""
    at_infinity(z::PVector{T, Int}, maxnorm)

Check whether `z` represents a point at infinity.
We declare `z` at infinity if the infinity norm of
its affine vector is larger than `maxnorm`.
"""
function at_infinity(z::PVector{<:Complex, Int}, maxnorm)
    at_inf = false
    @inbounds tol = maxnorm * maxnorm * abs2(z.data[z.homvar])
    @inbounds for k in 1:length(z)
        if k == z.homvar
            continue
        # we avoid the sqrt here by using the squared comparison
        elseif abs2(z.data[k]) > tol
            return true
        end
    end
    false
end
function at_infinity(z::PVector{<:Real, Int}, maxnorm)
    at_inf = false
    @inbounds tol = maxnorm * abs(z.data[z.homvar])
    @inbounds for k in 1:length(z)
        if k == z.homvar
            continue
        elseif abs(z.data[k]) > tol
            return true
        end
    end
    false
end

LinearAlgebra.norm(v::PVector, p::Real=2) = norm(v.data, p)
function LinearAlgebra.normalize!(v::PVector, p::Real=2)
    LinearAlgebra.normalize!(v.data, p)
    v
end

LinearAlgebra.dot(v::PVector, w::PVector) = LinearAlgebra.dot(v.data, w.data)

const VecView{T} = SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}
#
# """
#     ProdPVector(pvectors)
#
# Construct a product of `PVector`s. This
# """
# struct ProdPVector{T, H<:Union{Nothing, Int}} <: AbstractProjectiveVector{T}
#     data::Vector{T}
#     # It is faster to use a tuple instead but this puts
#     # additional stress on the compiler and makes things
#     # not inferable, so we do this not for now.
#     # But we should be able to change this easily later
#     pvectors::Vector{PVector{T, H, VecView{T}}}
# end
#
# function ProdPVector(vs)
#     data = copy(vs[1].data);
#     for k=2:length(vs)
#         append!(data, vs[k].data)
#     end
#     pvectors = _create_pvectors(data, vs)
#     ProdPVector(data, pvectors)
# end
#
# function (==)(v::ProdPVector, w::ProdPVector)
#     for (vᵢ, wᵢ) in zip(pvectors(v), pvectors(w))
#         if vᵢ != wᵢ
#             return false
#         end
#     end
#     true
# end
#
# function _create_pvectors(data, vs)
#     pvectors = Vector{PVector{eltype(data), VecView{eltype(data)}}}()
#     k = 1
#     for i = 1:length(vs)
#         n = length(vs[i])
#         vdata = view(data, k:k+n-1)
#         push!(pvectors, PVector(vdata, homvar(vs[i])))
#         k += n
#     end
#     pvectors
# end
#
# function Base.copy(z::ProdPVector)
#     data = copy(z.data)
#     new_pvectors = _create_pvectors(data, pvectors(z))
#     ProdPVector(data, new_pvectors)
# end
#
# function Base.similar(v::ProdPVector, ::Type{T}) where T
#     data = convert.(T, v.data)
#     ProdPVector(data, _create_pvectors(data, pvectors(v)))
# end
#
# """
#     pvectors(z::ProdPVector)
#
# Return the `PVector`s out of which the product `z` exists.
# """
# pvectors(z::ProdPVector) = z.pvectors
#
# """
#     raw(z::ProdPVector)
#
# access_input the underlying vector of the product `z`. Note that this
# is only a single vector. This is useful to pass
# the vector into some function which does not know the
# projective structure.
# """
# raw(v::ProdPVector) = v.data
#
# """
#     at_infinity(z::ProdPVector, maxnorm)
#
# Returns `true` if any vector of the product is at infinity.
# """
# function at_infinity(v::ProdPVector, maxnorm)
#     for vᵢ in pvectors(v)
#         if at_infinity(vᵢ, maxnorm)
#             return true
#         end
#     end
#     false
# end
#
# """
#     affine(z::ProdPVector)
#
# For each projective vector of the product return associated affine patch.
# """
# affine(z::ProdPVector) = affine.(z.pvectors)
#
# """
#     affine!(z::ProdPVector)
#
# Bring each projective vector of the product on the associated affine patch.
# """
# function affine!(v::ProdPVector)
#     for w in pvectors(v)
#         affine!(w)
#     end
#     v
# end
#
# function LinearAlgebra.normalize!(v::ProdPVector, p::Real=2)
#     for w in pvectors(v)
#         normalize!(w, p)
#     end
#     v
# end
#
# infinity_norm(z::ProdPVector) = maximum(infinity_norm, pvectors(z))
# function infinity_norm(v::ProdPVector, w::ProdPVector)
#     p₁ = pvectors(v)
#     p₂ = pvectors(w)
#     m = infinity_norm(p₁[1], p₂[1])
#     for k=2:length(p₁)
#         m = max(m, infinity_norm(p₁[k], p₂[k]))
#     end
#     m
# end


# AbstractVector interface

Base.size(z::AbstractProjectiveVector) = size(z.data)
Base.length(z::AbstractProjectiveVector) = length(z.data)
Base.getindex(z::AbstractProjectiveVector, i::Integer) = getindex(z.data, i)
Base.setindex!(z::AbstractProjectiveVector, zᵢ, i::Integer) = setindex!(z.data, zᵢ, i)
Base.lastindex(z::AbstractProjectiveVector) = lastindex(z.data)

end
