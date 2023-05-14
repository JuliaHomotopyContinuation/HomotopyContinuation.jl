export SystemHomotopy

"""
    SystemHomotopy(H::AbstactSystem)

Interprets system `F(x,p)` as the homotopy `H(x,t) = F(x, [t])`
"""
struct SystemHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    system::T
    t::Vector{ComplexF64}
    tt::TaylorVector{2,ComplexF64}
end
SystemHomotopy(F::AbstractSystem; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]) =
    SystemHomotopy(F, [ComplexF64(0)], TaylorVector{2}(ComplexF64, 1))

system_as_homotopy(F::AbstractSystem) = SystemHomotopy(F)

Base.size(H::SystemHomotopy) = size(H.system)

(H::SystemHomotopy)(x, t, p::Nothing = nothing) = H.system(x, [t])

on_affine_chart(H::SystemHomotopy, proj_dims = nothing;) =
    on_affine_chart(H.system, proj_dims)

function set_solution!(x::AbstractVector, H::SystemHomotopy, y::AbstractVector, t)
    set_solution!(x, H.system, y)
end

function ModelKit.evaluate!(u, H::SystemHomotopy, x, t, p::Nothing = nothing)
    H.t[1] = t
    evaluate!(u, H.system, x, H.t)
end
function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::SystemHomotopy,
    x,
    t,
    p::Nothing = nothing,
)
    H.t[1] = t
    evaluate_and_jacobian!(u, U, H.system, x, H.t)
end

function ModelKit.taylor!(u, v::Val, H::SystemHomotopy, tx, tt, p::Nothing = nothing)
    t, ṫ = tt
    H.tt[1, 1] = t
    H.tt[1, 2] = ṫ
    taylor!(u, v, H.system, tx, H.tt)
end

