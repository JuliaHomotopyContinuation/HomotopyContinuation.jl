export ModelKitSystem

"""
    ModelKitSystem(system::ModelKit.System) <: AbstractSystem
"""
struct ModelKitSystem{T} <: AbstractSystem
    system::ModelKit.CompiledSystem{T}
    degrees::Vector{Int}
end
#
function ModelKitSystem(system::ModelKit.System)
    csystem = ModelKit.compile(system)
    degrees = maxdegrees(system)
    ModelKitSystem(csystem, degrees)
end

Base.size(F::ModelKitSystem) = size(F.system)

cache(F::ModelKitSystem, x, p = nothing) = SystemNullCache()

@propagate_inbounds evaluate!(u, F::ModelKitSystem, x, ::SystemNullCache) =
    ModelKit.evaluate!(u, F.system, x)
evaluate(F::ModelKitSystem, x, ::SystemNullCache) = ModelKit.evaluate(F.system, x)
(F::ModelKitSystem)(x...) = evaluate(F, x...)

@propagate_inbounds jacobian!(U, F::ModelKitSystem, x, ::SystemNullCache) =
    ModelKit.jacobian!(U, F.system, x)
jacobian(F::ModelKitSystem, x, ::SystemNullCache) = ModelKit.jacobian(F.system, x)

@propagate_inbounds evaluate_and_jacobian!(u, U, F::ModelKitSystem, x, ::SystemNullCache) =
    ModelKit.evaluate_and_jacobian!(u, U, F.system, x)

degrees(F::ModelKitSystem) = F.degrees


# # parameter version
# @propagate_inbounds evaluate!(u, F::ModelKitSystem, x, p, ::SystemNullCache) =
#     ModelKit.evaluate!(u, F.system, x, p)
# evaluate(F::ModelKitSystem, x, p, ::SystemNullCache) = ModelKit.evaluate(F.system, x, p)
# @propagate_inbounds jacobian!(U, F::ModelKitSystem, x, p, ::SystemNullCache) =
#     ModelKit.jacobian!(U, F.system, x, p)
# jacobian(F::ModelKitSystem, x, p, ::SystemNullCache) = ModelKit.jacobian(F.system, x, p)
# @propagate_inbounds evaluate_jacobian!(u, U, F::ModelKitSystem, x, p, ::SystemNullCache) =
#     ModelKit.evaluate_jacobian!(u, U, F.system, x, p)
# evaluate_jacobian(F::ModelKitSystem, x, p, ::SystemNullCache) =
#     ModelKit.evaluate_jacobian(F.system, x, p)
