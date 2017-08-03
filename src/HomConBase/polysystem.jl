
# Constructors

PolySystem{T} = Vector{Poly{T}}

evaluate(F::PolySystem{T}, x::Vector{T}) where {T<:Number}= map(f -> evaluate(f, x), F)

@inline function jacobian(F::PolySystem{T}) where {T<:Number}
    DF = permutedims(hcat(gradient.(F)...), [2, 1])
    (x::Vector{T}) -> map(f -> evaluate(f, x), DF)
end

@inline is_homogenous(F::PolySystem) = all(is_homogenous, F)
@inline homogenize(F::PolySystem) = homogenize.(F)
@inline degrees(F::PolySystem) = deg.(F)
@inline nequations(F::PolySystem) = size(F, 1)
nvars(F::PolySystem) = nvars(F[1])

@inline function weyl_dot(F::PolySystem{T},G::PolySystem{T}) where {T<:Complex}
    sum(tup -> weyl_dot(tup[1],tup[2]), zip(F, G))
end

@inline function weyl_norm(F::PolySystem{T}) where {T<:Complex}
    √(weyl_dot(F, F))
end