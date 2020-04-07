"""
    ModelKitHomotopy(H::ModelKit.Homotopy, params = [])

Construct a `ModelKitHomotopy` with the given parameters fixed.
"""
struct ModelKitHomotopy{S,P<:Union{Nothing,AbstractVector}} <: AbstractHomotopy
    homotopy::ModelKit.CompiledHomotopy{S}
    parameters::P
end

ModelKitHomotopy(H::ModelKit.Homotopy, parameters = nothing) =
    ModelKitHomotopy(ModelKit.compile(H), parameters)
ModelKitHomotopy(H::ModelKit.CompiledHomotopy) = ModelKitHomotopy(H, nothing)

Base.size(H::ModelKitHomotopy) = size(H.homotopy)

evaluate!(u, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate!(u, H.homotopy, x, t, H.parameters)
evaluate_and_jacobian!(u, U, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.parameters)
function taylor!(u, v::Val, H::ModelKitHomotopy, x, t)
    ModelKit.taylor!(u, v, H.homotopy, x, t, H.parameters)
    u
end
