export SystemHomotopyComposition

"""
    SystemHomotopyComposition(F::AbstactSystem, H::AbstractHomotopy)

The homotopy ``J(x,t) = F(H(x,t))``.
"""
struct SystemHomotopyComposition{S<:AbstractSystem,T<:AbstractHomotopy} <: AbstractHomotopy
    F::S
    H::T
    composed::SystemHomotopy{CompositionSystem{HomotopySystem{T},S}}
end
SystemHomotopyComposition(
    F::AbstractSystem,
    H::AbstractHomotopy;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = SystemHomotopyComposition(F, H, system_as_homotopy(compose(F, homotopy_as_system(H))))

compose(F::AbstractSystem, H::AbstractHomotopy) = SystemHomotopyComposition(F, H)

Base.size(H::SystemHomotopyComposition) = size(H.composed)

(H::SystemHomotopyComposition)(x, t, p::Nothing = nothing) = H.composed(x, t)

function ModelKit.evaluate!(u, H::SystemHomotopyComposition, x, t, p::Nothing = nothing)
    evaluate!(u, H.composed, x, t)
end
function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::SystemHomotopyComposition,
    x,
    t,
    p::Nothing = nothing,
)
    evaluate_and_jacobian!(u, U, H.composed, x, t)
end

function ModelKit.taylor!(
    u,
    v::Val,
    H::SystemHomotopyComposition,
    tx,
    tt,
    p::Nothing = nothing,
)
    taylor!(u, v, H.composed, tx, tt)
end

