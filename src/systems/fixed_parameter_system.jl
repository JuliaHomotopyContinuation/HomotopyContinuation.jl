export FixedParameterSystem, fix_parameters

"""
    FixedParameterSystem(F:AbstractSystem, parameters)

Construct a system from the given [`AbstractSystem`](@ref) `F` with the given `parameters`
fixed.
"""
struct FixedParameterSystem{S<:AbstractSystem,T} <: AbstractSystem
    system::S
    parameters::Vector{T}
end
FixedParameterSystem(F::System, p) = FixedParameterSystem(ModelKitSystem(F), p)
Base.size(F::FixedParameterSystem) = size(F.system)

ModelKit.variables(F::FixedParameterSystem) = variables(F.system)
ModelKit.parameters(F::FixedParameterSystem) = Variable[]
ModelKit.variable_groups(F::FixedParameterSystem) = variable_groups(F.system)

(F::FixedParameterSystem)(x, p = nothing) = F.system(x, F.parameters)

evaluate!(u, F::FixedParameterSystem, x, p = nothing) =
    evaluate!(u, F.system, x, F.parameters)
evaluate_and_jacobian!(u, U, F::FixedParameterSystem, x, p = nothing) =
    evaluate_and_jacobian!(u, U, F.system, x, F.parameters)
taylor!(u, v::Val, F::FixedParameterSystem, tx, p = nothing) =
    taylor!(u, v, F.system, tx, F.parameters)

"""
    fix_parameters(F::AbstractSystem, p)

Fix the parameters of the given system `F`. Returns a [`FixedParameterSystem`](@ref).
"""
fix_parameters(F::AbstractSystem, p) = FixedParameterSystem(F, p)
