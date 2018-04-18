module ProjectiveVectors

import Base: ==
import ..Utilities: infinity_norm

export AbstractProjectiveVector,
    PVector,
    ProdPVector,
    raw,
    homvar,
    affine,
    affine!,
    embed,
    at_infinity,
    raw,
    pvectors,
    infinity_norm,
    infinity_norm_fast,
    unsafe_infinity_norm,
    converteltype

abstract type AbstractProjectiveVector{T} <: AbstractVector{T} end

"""
    PVector(z, homvar)

A `PVector` represents a projective vector `z` which is the result
of the embedding of an affine vector `x` into projective space.
Note that the constructor *does not* perform the embedding but rather
is a simple wrapper. In order to embed an affine vector use
[`embed`](@ref).
"""
struct PVector{T, V<:AbstractVector{T}} <: AbstractProjectiveVector{T}
    data::V
    homvar::Int

    function PVector{T, V}(data, homvar) where {T, V}
        @assert 1 ≤ homvar ≤ length(data)
        new(data, homvar)
    end
end
PVector(x::V, homvar::Int) where {T, V<:AbstractVector{T}} = PVector{T, V}(x, homvar)

Base.copy(z::PVector{T, V}) where {T, V} = PVector{T, V}(copy(z.data), z.homvar)

"""
    raw(z::PVector)

Access the underlying vector of `z`. This is useful to pass
the vector into some function which does not know the
projective structure.
"""
raw(z::PVector) = z.data

"""
    homvar(z::PVector)

Get the index of the homogenous variable.
"""
homvar(z::PVector) = z.homvar

(==)(v::PVector, w::PVector) = v.homvar == w.homvar && v.data == w.data

Base.copy(z::PVector) = PVector(copy(z.data), z.homvar)

"""
    converteltype(v::PVector, T)
"""
converteltype(v::PVector, ::Type{T}) where T = PVector(convert.(T, v.data), v.homvar)
converteltype(v::PVector{T}, ::Type{T}) where T = copy(v)

"""
    embed(x::Vector, homvar)::PVector
"""
function embed(x::Vector{T}, homvar) where T
    k = 1
    data = Vector{T}(length(x) + 1)
    for k in 1:length(data)
        if k == homvar
            data[k] = one(T)
        else
            i = k < homvar ? k : k - 1
            data[k] = x[i]
        end
    end
    PVector(data, homvar)
end

"""
    affine_normalizer(z::PVector)

Returns a scalar to bring `z` on its affine patch.
"""
affine_normalizer(z::PVector) = inv(z.data[z.homvar])
"""
    abs2_ffine_normalizer(z::PVector)

Returns the squared absolute value of the normalizer.
"""
abs2_affine_normalizer(z::PVector) = inv(abs2(z.data[z.homvar]))

"""
    affine(z::PVector)::Vector

Return the corresponding affine vector.
"""
function affine(z::PVector{T}) where {T}
    x = Vector{T}(length(z) - 1)
    normalizer = affine_normalizer(z)
    for k in 1:length(z)
        if k == z.homvar
            continue
        end
        i = k < z.homvar ? k : k - 1
        x[i] = z.data[k] * normalizer
    end
    x
end

"""
    affine!(z::PVector)

Bring the projective vector on associated affine patch.
"""
affine!(z::PVector) = scale!(z.data, affine_normalizer(z))

"""
    infinity_norm(z::PVector)

Compute the ∞-norm of `z`. If `z` is a complex vector this is more efficient
than `norm(z, Inf)`.

    infinity_norm(z₁::PVector, z₂::PVector)

Compute the ∞-norm of `z₁-z₂` by bringing both vectors first on their respective
affine patch. This therefore only makes sense if both vectors have the
same affine patch.
"""
function infinity_norm(z::PVector{<:Complex})
    sqrt(maximum(abs2, raw(z)) * abs2_affine_normalizer(z))
end

function infinity_norm(z₁::PVector{<:Complex}, z₂::PVector{<:Complex})
    normalizer₁ = affine_normalizer(z₁)
    normalizer₂ = -affine_normalizer(z₂)
    m = abs2(muladd(z₁[1], normalizer₁, z₂[1] * normalizer₂))
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
function unsafe_infinity_norm(z₁::PVector{<:Complex}, z₂::PVector{<:Complex})
    m = abs2(z₁[1] - z₂[1])
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
    at_infinity(z::PVector, maxnorm)

Check whether `z` represents a point at infinity.
We declare `z` at infinity if the infinity norm of
its affine vector is larger than `maxnorm`.
"""
function at_infinity(z::PVector{<:Complex}, maxnorm)
    at_inf = false
    tol = maxnorm * maxnorm * abs2(z.data[z.homvar])
    for k in 1:length(z)
        if k == z.homvar
            continue
        # we avoid the sqrt here by using the squared comparison
        elseif abs2(z.data[k]) > tol
            return true
        end
    end
    false
end
function at_infinity(z::PVector{<:Real}, maxnorm)
    at_inf = false
    tol = maxnorm * abs(z.data[z.homvar])
    for k in 1:length(z)
        if k == z.homvar
            continue
        elseif abs(z.data[k]) > tol
            return true
        end
    end
    false
end

Base.LinAlg.norm(v::PVector, p::Real=2) = norm(v.data, p)
function Base.LinAlg.normalize!(v::PVector, p::Real=2)
    normalize!(v.data, p)
    v
end


const VecView{T} = SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}

"""
    ProdPVector(pvectors)

Construct a product of `PVector`s. This
"""
struct ProdPVector{T} <: AbstractProjectiveVector{T}
    data::Vector{T}
    # It is faster to use a tuple instead but this puts
    # additional stress on the compiler and makes things
    # not inferable, so we do this not for now.
    # But we should be able to change this easily later
    pvectors::Vector{PVector{T, VecView{T}}}
end

function ProdPVector(vs)
    data = copy(vs[1].data);
    for k=2:length(vs)
        append!(data, vs[k].data)
    end
    pvectors = _create_pvectors(data, vs)
    ProdPVector(data, pvectors)
end

function (==)(v::ProdPVector, w::ProdPVector)
    for (vᵢ, wᵢ) in zip(pvectors(v), pvectors(w))
        if vᵢ != wᵢ
            return false
        end
    end
    true
end

function _create_pvectors(data, vs)
    pvectors = Vector{PVector{eltype(data), VecView{eltype(data)}}}()
    k = 1
    for i = 1:length(vs)
        n = length(vs[i])
        vdata = view(data, k:k+n-1)
        push!(pvectors, PVector(vdata, homvar(vs[i])))
        k += n
    end
    pvectors
end

function Base.copy(z::ProdPVector)
    data = copy(z.data)
    new_pvectors = _create_pvectors(data, pvectors(z))
    ProdPVector(data, new_pvectors)
end

"""
    converteltype(v::ProdPVector, T)
"""
function converteltype(v::ProdPVector, ::Type{T}) where T
    data = convert.(T, v.data)
    pvectors = _create_pvectors(data, vs)
    ProdPVector(data, pvectors)
end
converteltype(v::ProdPVector{T}, ::Type{T}) where T = v

"""
    pvectors(z::ProdPVector)

Return the `PVector`s out of which the product `z` exists.
"""
pvectors(z::ProdPVector) = z.pvectors

"""
    raw(z::ProdPVector)

Access the underlying vector of the product `z`. Note that this
is only a single vector. This is useful to pass
the vector into some function which does not know the
projective structure.
"""
raw(v::ProdPVector) = v.data

"""
    at_infinity(z::ProdPVector, maxnorm)

Returns `true` if any vector of the product is at infinity.
"""
function at_infinity(v::ProdPVector, maxnorm)
    for vᵢ in pvectors(v)
        if at_infinity(vᵢ, maxnorm)
            return true
        end
    end
    false
end

"""
    affine(z::ProdPVector)

For each projective vector of the product return associated affine patch.
"""
affine(z::ProdPVector) = affine.(v.pvectors)

"""
    affine!(z::ProdPVector)

Bring each projective vector of the product on the associated affine patch.
"""
function affine!(v::ProdPVector)
    for w in pvectors(v)
        affine!(w)
    end
    v
end

function Base.LinAlg.normalize!(v::ProdPVector, p::Real=2)
    for w in pvectors(v)
        normalize!(w, p)
    end
    v
end

infinity_norm(z::ProdPVector) = maximum(infinity_norm, pvectors(z))
function infinity_norm(v::ProdPVector, w::ProdPVector)
    p₁ = pvectors(v)
    p₂ = pvectors(w)
    m = infinity_norm(p₁[1], p₂[1])
    for k=2:length(p1)
        m = max(m, infinity_norm(p₁[k], p₂[k]))
    end
    m
end


# AbstractVector interface

Base.size(z::AbstractProjectiveVector) = size(z.data)
Base.length(z::AbstractProjectiveVector) = length(z.data)
Base.getindex(z::AbstractProjectiveVector, i::Integer) = getindex(z.data, i)
Base.setindex!(z::AbstractProjectiveVector, zᵢ, i::Integer) = setindex!(z.data, zᵢ, i)
Base.endof(z::AbstractProjectiveVector) = endof(z.data)
Base.eltype(z::AbstractProjectiveVector{T}) where T = T

end
