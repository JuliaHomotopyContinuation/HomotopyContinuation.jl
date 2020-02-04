"""
    ModelKitHomotopy(H::ModelKit.Homotopy, params = [])

Construct a `ModelKitHomotopy` with the given parameters fixed.
"""
struct ModelKitHomotopy{S,T} <: AbstractHomotopy
    homotopy::ModelKit.CompiledHomotopy{S}
    parameters::Vector{T}
end

ModelKitHomotopy(H::ModelKit.Homotopy, parameters = ComplexF64[]) =
    ModelKitHomotopy(ModelKit.compile(H), parameters)
ModelKitHomotopy(H::ModelKit.CompiledHomotopy) =
    ModelKitHomotopy(H, ComplexF64[])

Base.size(H::ModelKitHomotopy) = size(H.homotopy)

evaluate!(u, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate!(u, H.homotopy, x, t, H.parameters)
# jacobian!(U, H::ModelKitHomotopy, x, t) =
#     ModelKit.jacobian!(U, H.homotopy, x, t, H.parameters)
evaluate_and_jacobian!(u, U, H::ModelKitHomotopy, x, t) =
    ModelKit.evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.parameters)
diff_t!(u, H::ModelKitHomotopy, x, t, dx = ()) =
    ModelKit.diff_t!(u, H.homotopy, x, t, dx, H.parameters)
