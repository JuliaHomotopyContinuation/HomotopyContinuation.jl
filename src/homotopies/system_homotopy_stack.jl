export SystemHomotopyStack

"""
    SystemHomotopyStack(F::Union{AbstractSystem,System}; start_parameters, target_parameters)
    SystemHomotopyStack(F::Union{AbstractSystem,System}, start_parameters, target_parameters)

Construct the parameter homotopy ``H(x,t) = F(x; t p + (1 - t) q)`` where ``p`` is
`start_parameters` and ``q`` is `target_parameters`.
"""
struct SystemHomotopyStack{S<:AbstractSystem,T<:AbstractHomotopy} <: AbstractHomotopy
    F::S
    H::T
end

stack(F::AbstractSystem, H::AbstractHomotopy) = SystemHomotopyStack(F, H)
homotopy(H::SystemHomotopyStack) = H.H
Base.size(C::SystemHomotopyStack) = (size(C.F, 1) + size(C.H, 1), size(C.F, 2))

function ModelKit.evaluate!(
    u,
    C::SystemHomotopyStack,
    x::AbstractVector,
    t,
    p::Nothing = nothing,
)
    m_f = size(C.F, 1)
    u_f = @view u[1:m_f]
    evaluate!(u_f, C.F, x, p)

    m_h = size(C.H, 1)
    u_h = @view u[(m_f+1):(m_f+m_h)]
    evaluate!(u_h, C.H, x, t, p)
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    C::SystemHomotopyStack,
    x,
    t,
    p::Nothing = nothing,
)
    m_f, n = size(C.F)
    u_f = @view u[1:m_f]
    U_f = @view U[1:m_f, 1:n]
    evaluate_and_jacobian!(u_f, U_f, C.F, x, p)

    m_h = size(C.H, 1)
    u_h = @view u[(m_f+1):(m_f+m_h)]
    U_h = @view U[(m_f+1):(m_f+m_h), 1:n]
    evaluate_and_jacobian!(u_h, U_h, C.H, x, t, p)

    nothing
end

function ModelKit.taylor!(u, v::Val, C::SystemHomotopyStack, tx, t, p::Nothing = nothing)
    m_f = size(C.F, 1)
    u_f = @view u[1:m_f]
    taylor!(u_f, v, C.F, tx, p)

    m_h = size(C.H, 1)
    u_h = @view u[(m_f+1):(m_f+m_h)]
    taylor!(u_h, v, C.H, tx, t, p)

    u
end

