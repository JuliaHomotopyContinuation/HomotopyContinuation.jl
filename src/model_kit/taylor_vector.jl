"""
    TaylorVector{N,T}

A data structure representing a vector of Taylor series with `N` terms.
Each element is an `NTuple{N,T}`.

    TaylorVector(T, N::Integer, n::Integer)

Create a vector of `n` `NTuple{N,T}`s.
"""
struct TaylorVector{N,T} <: AbstractVector{NTuple{N,T}}
    data::LinearAlgebra.Transpose{T,Matrix{T}}
    views::NTuple{N,SubArray{T,1,Matrix{T},Tuple{Int,UnitRange{Int}},true}}
end
TaylorVector{N}(TV::TaylorVector{M,T}) where {N,M,T} =
    TaylorVector{N,T}(TV.data, TV.views[1:N])
function TaylorVector{N}(data::Matrix) where {N}
    views = tuple((view(data, i, 1:size(data, 2)) for i = 1:N)...)
    TaylorVector(LinearAlgebra.transpose(data), views)
end
function TaylorVector{N}(T, n::Integer) where {N}
    TaylorVector{N}(zeros(T, N, n))
end

function Base.show(io::IO, ::MIME"text/plain", X::TaylorVector)
    summary(io, X)
    isempty(X) && return
    println(io, ":")
    Base.print_array(io, map(i -> X[i], 1:length(X)))
end

Base.length(TV::TaylorVector) = size(TV.data, 1)
Base.size(TV::TaylorVector) = (length(TV),)
Base.eltype(::Type{TaylorVector{N,T}}) where {N,T} = NTuple{N,T}
Base.IndexStyle(::TaylorVector) = Base.IndexLinear()

"""
    vectors(TV::TaylorVec{N})

Return the Taylor series as `N` seperate vectors.
"""
vectors(TV::TaylorVector) = TV.views

@generated function Base.getindex(TV::TaylorVector{N}, i::Integer) where {N}
    quote
        Base.@_propagate_inbounds_meta
        x = TV.data
        $(Expr(:tuple, (:(x[i, $k]) for k = 1:N)...))
    end
end
Base.getindex(TV::TaylorVector, i::Integer, j::Integer) = getindex(TV.data, i, j)

function Base.setindex!(TV::TaylorVector{N,T}, x, i::Integer) where {N,T}
    setindex!(TV, convert(NTuple{N,T}, x), i)
end
function Base.setindex!(TV::TaylorVector{1}, x::Number, i::Integer)
    TV.data[i] = x
    x
end
@generated function Base.setindex!(
    TV::TaylorVector{N,T},
    x::NTuple{N,T},
    i::Integer,
) where {N,T}
    quote
        Base.@_propagate_inbounds_meta
        d = TV.data
        $(Expr(:tuple, (:(d[i, $k]) for k = 1:N)...)) = x
        x
    end
end
Base.setindex!(TV::TaylorVector, x, i::Integer, j::Integer) = setindex!(TV.data, x, i, j)
