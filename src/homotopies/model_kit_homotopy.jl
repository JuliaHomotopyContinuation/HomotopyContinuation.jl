export ModelKitHomotopy

"""
    ModelKitHomotopy(H::Homotopy, parameters = nothing)

Construct a homotopy from the given homotopy `H` with the given `parameters` fixed.
The difference to `H` is that this compiles a straight line programm for the fast
evaluation of `H` and that `ModelKitHomotopy <: AbstractHomotopy`.
"""
struct ModelKitHomotopy{S,P<:Union{Nothing,AbstractVector}} <: AbstractHomotopy
    homotopy::ModelKit.CompiledHomotopy{S}
    parameters::P
end

ModelKitHomotopy(H::ModelKit.Homotopy, parameters = nothing) =
    ModelKitHomotopy(ModelKit.compile(H), parameters)
ModelKitHomotopy(H::ModelKit.CompiledHomotopy, parameters = nothing) =
    ModelKitHomotopy(H, parameters)

Base.size(H::ModelKitHomotopy) = size(H.homotopy)

evaluate!(u, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate!(u, H.homotopy, x, t, H.parameters)
evaluate_and_jacobian!(u, U, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.parameters)
function taylor!(u, v::Val, H::ModelKitHomotopy, x, t)
    ModelKit.taylor!(u, v, H.homotopy, x, t, H.parameters)
    u
end
