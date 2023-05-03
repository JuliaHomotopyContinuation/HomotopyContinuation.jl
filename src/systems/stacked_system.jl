export StackedSystem, StackedSystem

"""
    StackedSystem(F::AbstractSystem, G::AbstractSystem)

Construct the system ``[F(x;p); G(x;p)]``.
"""
struct StackedSystem{S1<:AbstractSystem,S2<:AbstractSystem} <: AbstractSystem
    f::S1
    g::S2
    parameters::Vector{Variable}
end
function StackedSystem(g::AbstractSystem, f::AbstractSystem)
    variables(g) == variables(f) || throw(
        ArgumentError(
            "Cannot create the stacked system [F;g] F since the variabels of `F` and `G` don't match.",
        ),
    )

    pg = parameters(g)
    pf = parameters(f)
    (isempty(pf) || isempty(pg) || pf == pg) || throw(
        ArgumentError(
            "Cannot construct a stack of two system with different sets of parameters.",
        ),
    )
    p = isempty(pg) ? [pg; pf] : pg

    StackedSystem(g, f, p)
end
StackedSystem(
    g::AbstractSystem,
    f::System;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = StackedSystem(g, fixed(f; compile = compile))
StackedSystem(
    g::System,
    f::AbstractSystem;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = StackedSystem(fixed(g; compile = compile), f)
StackedSystem(g::System, f::System; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]) =
    StackedSystem(fixed(g; compile = compile), fixed(f; compile = compile))

stack(g, f) = StackedSystem(g, f)

Base.size(C::StackedSystem) = (size(C.f, 1) + size(C.g, 1), size(C.f, 2))

ModelKit.variables(F::StackedSystem) = variables(F.f)
ModelKit.parameters(F::StackedSystem) = F.parameters
ModelKit.variable_groups(F::StackedSystem) = variable_groups(F.f)

function Base.show(io::IO, C::StackedSystem)
    println(io, "Composition G ∘ F:")
    println(io, "F: ")
    show(io, C.f)
    println(io, "\nG: ")
    show(io, C.g)
end

(C::StackedSystem)(x, p = nothing) = [C.f(x, p); C.g(x, p)]
function ModelKit.evaluate!(u, C::StackedSystem, x::AbstractVector, p = nothing)
    m_f = size(C.f, 1)
    u_f = @view u[1:m_f]
    evaluate!(u_f, C.f, x, p)

    m_g = size(C.g, 1)
    u_g = @view u[(m_f+1):(m_f+m_g)]
    evaluate!(u_g, C.g, x, p)
    u
end

function ModelKit.evaluate_and_jacobian!(u, U, C::StackedSystem, x, p = nothing)
    m_f, n = size(C.f)
    u_f = @view u[1:m_f]
    U_f = @view U[1:m_f, 1:n]
    evaluate_and_jacobian!(u_f, U_f, C.f, x, p)

    m_g = size(C.g, 1)
    u_g = @view u[(m_f+1):(m_f+m_g)]
    U_g = @view U[(m_f+1):(m_f+m_g), 1:n]
    evaluate_and_jacobian!(u_g, U_g, C.g, x, p)

    nothing
end

function ModelKit.taylor!(u, v::Val, C::StackedSystem, tx, p = nothing)
    m_f = size(C.f, 1)
    u_f = @view u[1:m_f]
    taylor!(u_f, v, C.f, tx, p)

    m_g = size(C.g, 1)
    u_g = @view u[(m_f+1):(m_f+m_g)]
    taylor!(u_g, v, C.g, tx, p)

    u
end
