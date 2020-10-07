"""
    TaylorVector{N,T}

A data structure representing a vector of Taylor series with `N` terms.
Each element is an `NTuple{N,T}`.

    TaylorVector(T, N::Integer, n::Integer)

Create a vector of `n` `NTuple{N,T}`s.
"""
struct TaylorVector{N,T} <: AbstractVector{NTuple{N,T}}
    data::Matrix{T}
end
TaylorVector{N}(data::Matrix{T}) where {N,T} = TaylorVector{N,T}(data)
TaylorVector{N}(TV::TaylorVector{M,T}) where {N,M,T} = TaylorVector{N,T}(TV.data)
TaylorVector{N}(T, n::Integer) where {N} = TaylorVector{N}(zeros(T, N, n))

function Base.show(io::IO, ::MIME"text/plain", X::TaylorVector)
    summary(io, X)
    isempty(X) && return
    println(io, ":")
    Base.print_array(io, map(i -> X[i], 1:length(X)))
end

Base.length(TV::TaylorVector) = size(TV.data, 2)
Base.size(TV::TaylorVector) = (length(TV),)
Base.eltype(::Type{TaylorVector{N,T}}) where {N,T} = NTuple{N,T}
Base.IndexStyle(::TaylorVector) = Base.IndexLinear()

Base.fill!(v::TaylorVector, x) = fill!(v.data, x)

"""
    vectors(TV::TaylorVec{N})

Return the Taylor series as `N` seperate vectors.
"""
@generated function vectors(TV::TaylorVector{N}) where {N}
    Expr(:tuple, (:(view(TV.data, $k, :)) for k = 1:N)...)
end

@generated function Base.getindex(TV::TaylorVector{N}, i::Integer) where {N}
    quote
        Base.@_propagate_inbounds_meta
        x = TV.data
        $(Expr(:tuple, (:(x[$k, i]) for k = 1:N)...))
    end
end
Base.getindex(TV::TaylorVector, i::Integer, j::Integer) = getindex(TV.data, j, i)

@generated function Base.setindex!(
    TV::TaylorVector{N,T},
    x::NTuple{M,S},
    i::Integer,
) where {N,T,M,S}
    tup = Expr(:tuple, (i ≤ M ? :(convert(T, x[$i])) : :(zero(T)) for i = 1:N)...)
    :(setindex!(TV, $tup, i))
end
Base.setindex!(TV::TaylorVector, x, i::Integer) = setindex!(TV, make_ntuple(x), i)
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
        $(Expr(:tuple, (:(d[$k, i]) for k = 1:N)...)) = x
        x
    end
end
Base.setindex!(TV::TaylorVector, x, i::Integer, j::Integer) = setindex!(TV.data, x, j, i)
