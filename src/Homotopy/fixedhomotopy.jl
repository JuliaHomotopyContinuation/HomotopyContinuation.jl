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

homogenize(H::FixedHomotopy) = FixedHomotopy(homogenize(H.homotopy), H.t)
homogenized(H::FixedHomotopy) = homogenized(H.homotopy)
ishomogenous(H::FixedHomotopy) = ishomogenous(H.homotopy)
nvariables(H::FixedHomotopy) = nvariables(H.homotopy)
degrees(H::FixedHomotopy) = degrees(H.homotopy)

Base.length(H::FixedHomotopy) = length(H.homotopy)
Base.eltype(::FixedHomotopy{T, H, S}) where {T, H, S} = T
