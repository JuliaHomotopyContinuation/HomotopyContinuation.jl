export ModelKitSystem

"""
    ModelKitSystem(F:System, parameters = nothing)

Construct a system from the given [`System`](@ref) `F` with the given `parameters` fixed.
The difference to `F` is that this compiles a straight line programm for the fast
evaluation of `F` and that `ModelKitSystem <: AbstractSystem`.
"""
struct ModelKitSystem{S,T<:Union{Nothing,AbstractVector}} <: AbstractSystem
    system::ModelKit.CompiledSystem{S}
    parameters::T
end

ModelKitSystem(F::ModelKit.CompiledSystem) = ModelKitSystem(F, nothing)
function ModelKitSystem(F::ModelKit.System, p = nothing)
    if p !== nothing
        length(p) === nparameters(F) || throw(ArgumentError("Length of provided parameters doesn't match the number of parameters of the system."))
    end
    ModelKitSystem(ModelKit.compile(F), p)
end

Base.size(F::ModelKitSystem) = size(F.system)

function Base.show(io::IO, F::ModelKitSystem)
    println(io, typeof(F), ":")
    println(io, F.system)
    if F.parameters !== nothing
        println(io, "Parameters:")
        show(io, F.parameters)
    end
end

(F::ModelKitSystem{T,Nothing})(x, p = nothing) where {T} = F.system.system(x, p)
(F::ModelKitSystem{T,<:AbstractVector})(x) where {T} = F.system.system(x, F.parameters)
evaluate!(u, F::ModelKitSystem{T,<:AbstractVector}, x, p::Nothing = nothing) where {T} =
    ModelKit.evaluate!(u, F.system, x, F.parameters)
evaluate!(u, F::ModelKitSystem{T,Nothing}, x, p = nothing) where {T} =
    ModelKit.evaluate!(u, F.system, x, p)

evaluate_and_jacobian!(
    u,
    U,
    F::ModelKitSystem{T,<:AbstractVector},
    x,
    p::Nothing = nothing,
) where {T} = ModelKit.evaluate_and_jacobian!(u, U, F.system, x, F.parameters)
evaluate_and_jacobian!(u, U, F::ModelKitSystem{T,Nothing}, x, p = nothing) where {T} =
    ModelKit.evaluate_and_jacobian!(u, U, F.system, x, p)

taylor!(
    u,
    v::Val,
    F::ModelKitSystem{T,<:AbstractVector},
    tx::TaylorVector,
    p::Nothing = nothing,
) where {T} = ModelKit.taylor!(u, v, F.system, tx, F.parameters)
taylor!(u, v::Val, F::ModelKitSystem{T,Nothing}, tx::TaylorVector, p = nothing) where {T} =
    ModelKit.taylor!(u, v, F.system, tx, p)
