"""
    FixedHomotopy(homotopy, t)

Fix the homotopy at `t`, i.e. consider ``H(x,t)`` as ``H(x)``
"""
struct FixedHomotopy{T<:Number, H<:AbstractHomotopy{T}, S<:Number} <: AbstractPolySystem{T}
    homotopy::H
    t::S
end
function FixedHomotopy(hom::H, t::S) where {T<:Number, H<:AbstractHomotopy{T}, S<:Number}
    FixedHomotopy{T, H, S}(hom, t)
end

@inline evaluate(H::FixedHomotopy, x) = evaluate(H.homotopy, x, H.t)
(H::FixedHomotopy)(x) = evaluate(H,x)

function Base.promote_rule(::Type{FixedHomotopy{T1,H1,S1}}, ::Type{FixedHomotopy{T2, H2, S2}}) where {T1, T2, H1, H2, S1, S2}
    FixedHomotopy{promote_type(T1,T2), promote_type(H1,H2), promote_type(S1,S2)}
end

function Base.promote_rule(::Type{FixedHomotopy{T1,H1,S1}}, ::Type{T2}) where {T1, T2, H1, S1}
    FixedHomotopy{promote_type(T1,T2), promote_type(H1,T2), promote_type(S1,T2)}
end

function Base.convert(::Type{FixedHomotopy{T1,H1,S1}}, hom::FixedHomotopy) where{T1, H1, S1}
    FixedHomotopy(convert(H1, hom.homotopy), convert(S1, hom.t))
end

function Base.show(io::IO, H::FixedHomotopy)
    print(io, typeof(H), ":\n")
    println(io, "* homotopy:")
    println(io, H.homotopy)
    println(io, "* t: ", H.t)
end
function differentiate(H::FixedHomotopy)
    dh = differentiate(H.homotopy)
    (x) -> dh(x, H.t)
end

substitute(H::FixedHomotopy, pair) = FixedHomotopy(substitute(H.homotopy, pair), H.t)
homogenize(H::FixedHomotopy) = FixedHomotopy(homogenize(H.homotopy), H.t)
homogenized(H::FixedHomotopy) = homogenized(H.homotopy)
ishomogenous(H::FixedHomotopy) = ishomogenous(H.homotopy)
nvariables(H::FixedHomotopy) = nvariables(H.homotopy)
degrees(H::FixedHomotopy) = degrees(H.homotopy)

Base.length(H::FixedHomotopy) = length(H.homotopy)
coefftype(::FixedHomotopy{T, H, S}) where {T, H, S} = T
