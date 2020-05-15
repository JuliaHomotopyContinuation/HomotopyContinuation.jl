export ModelKitSystem

"""
    ModelKitSystem(F:System)

Construct a system from the given [`System`](@ref) `F`.
The difference to `F` is that this compiles a straight line programm for the fast
evaluation of `F` and that `ModelKitSystem <: AbstractSystem`.
"""
struct ModelKitSystem{S} <: AbstractSystem
    system::ModelKit.CompiledSystem{S}
end
ModelKitSystem(F::ModelKit.System) = ModelKitSystem(ModelKit.compile(F))

Base.size(F::ModelKitSystem) = size(F.system)
ModelKit.variables(F::ModelKitSystem) = variables(F.system.system)
ModelKit.parameters(F::ModelKitSystem) = parameters(F.system.system)
ModelKit.variable_groups(F::ModelKitSystem) = variable_groups(F.system.system)

function Base.show(io::IO, F::ModelKitSystem)
    println(io, typeof(F), ":")
    println(io, F.system)
end

(F::ModelKitSystem)(x, p = nothing) = F.system.system(x, p)
evaluate!(u, F::ModelKitSystem, x, p = nothing) = ModelKit.evaluate!(u, F.system, x, p)
evaluate_and_jacobian!(u, U, F::ModelKitSystem, x, p = nothing) =
    ModelKit.evaluate_and_jacobian!(u, U, F.system, x, p)
taylor!(u, v::Val, F::ModelKitSystem, tx::TaylorVector, p = nothing) =
    ModelKit.taylor!(u, v, F.system, tx, p)
