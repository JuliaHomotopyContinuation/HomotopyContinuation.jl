export HomotopySystem, fix_parameters

"""
    HomotopySystem(H::AbstractHomotopy)

Interprets a homotopy `H(x,t)` as a system `F(x,p) = H(x,p[1])` with a single parameter `p[1] = t`.
"""
struct HomotopySystem{H<:AbstractHomotopy} <: AbstractSystem
    homotopy::H

    HomotopySystem(
        H::T;
        compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    ) where {T<:AbstractHomotopy} = new{T}(H)
end


homotopy_as_system(H::AbstractHomotopy) = HomotopySystem(H)

Base.size(F::HomotopySystem) = size(F.homotopy)

ModelKit.variables(F::HomotopySystem) = variables(F.homotopy)
ModelKit.parameters(F::HomotopySystem) = [Variable(:t)]
ModelKit.variable_groups(F::HomotopySystem) = variable_groups(F.homotopy)

(F::HomotopySystem)(x, p::AbstractVector) = F.homotopy(x, p[1])

ModelKit.evaluate!(u, F::HomotopySystem, x, p::AbstractVector) =
    evaluate!(u, F.homotopy, x, p[1])
ModelKit.evaluate_and_jacobian!(u, U, F::HomotopySystem, x, p::AbstractVector) =
    evaluate_and_jacobian!(u, U, F.homotopy, x, p[1])
function ModelKit.taylor!(u, v::Val, F::HomotopySystem, tx, p::AbstractVector)
    taylor!(u, v, F.homotopy, tx, (p[1, 1], p[1, 2]))
end

