abstract type AbstractSystem end

struct ModelKitSystem{S,T<:Union{Nothing,AbstractVector}} <: AbstractSystem
    system::ModelKit.CompiledSystem{S}
    parameters::T
end

ModelKitSystem(F::ModelKit.CompiledSystem) = ModelKitSystem(F, nothing)
ModelKitSystem(F::ModelKit.System, p = nothing) = ModelKitSystem(ModelKit.compile(F), p)

Base.size(F::ModelKitSystem) = size(F.system)

evaluate!(u, F::ModelKitSystem, x) = ModelKit.evaluate!(u, F.system, x, F.parameters)
evaluate!(u, F::ModelKitSystem{<:Any,Nothing}, x, p = nothing) =
    ModelKit.evaluate!(u, F.system, x, p)

evaluate_and_jacobian!(u, U, F::ModelKitSystem, x) =
    ModelKit.evaluate_and_jacobian!(u, U, F.system, x, F.parameters)
evaluate_and_jacobian!(u, U, F::ModelKitSystem{<:Any,Nothing}, x, p = nothing) =
    ModelKit.evaluate_and_jacobian!(u, U, F.system, x, p)

diff_t!(u, v::Val, F::ModelKitSystem, x, dx = ()) =
    ModelKit.diff_t!(u, v, F.system, x, dx, F.parameters)
diff_t!(u, v::Val, F::ModelKitSystem{<:Any,Nothing}, x, dx = (), p = nothing, dp = ()) =
    ModelKit.diff_t!(u, v, F.system, x, dx, p, dp)

taylor!(u, v::Val, F::ModelKitSystem, x, dx = ()) =
    ModelKit.taylor!(u, v, F.system, x, dx, F.parameters)
taylor!(u, v::Val, F::ModelKitSystem{<:Any,Nothing}, x, dx = (), p = nothing, dp = ()) =
    ModelKit.taylor!(u, v, F.system, x, dx, p, dp)
