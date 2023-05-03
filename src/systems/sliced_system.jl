export SlicedSystem, slice

"""
    SlicedSystem(F::AbstractSystem, L::LinearSubspace)

Given a polynomial system `F` with the affine linear subspace `L`
defined by ``Ax = b`` this constructs the system `[F(x); A x - b] = 0`.
See also [`slice`](@ref).
"""
struct SlicedSystem{S<:AbstractSystem,T} <: AbstractSystem
    system::S
    affine::LinearSubspace{T}
    sliced::StackedSystem{S,LinearSystem}
end
SlicedSystem(
    F::System,
    A::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = SlicedSystem(fixed(F; compile = compile), A)
SlicedSystem(
    F::AbstractSystem,
    A::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = SlicedSystem(F, A, StackedSystem(F, LinearSystem(A; variables = variables(F))))

"""
    slice(F::System, L::LinearSubspace; compile = $(COMPILE_DEFAULT[]))
    slice(F::AbstractSystem, L::LinearSubspace)

Slice the zero set of the polynomial system `F` with the affine linear subspace `L`
defined by ``Ax = b``. This constructs the system `[F(x); A x - b] = 0`.
See also [`SlicedSystem`](@ref).
"""
slice(F::System, A::LinearSubspace; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]) =
    SlicedSystem(fixed(F; compile = compile), A)
slice(
    F::AbstractSystem,
    A::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = SlicedSystem(F, A)
Base.intersect(F::Union{System,AbstractSystem}, A::LinearSubspace) = SlicedSystem(F, A)
Base.intersect(A::LinearSubspace, F::Union{System,AbstractSystem}) = SlicedSystem(F, A)

ModelKit.variables(F::SlicedSystem) = variables(F.system)
ModelKit.parameters(F::SlicedSystem) = parameters(F.system)
ModelKit.variable_groups(F::SlicedSystem) = nothing

linear_subspace(A::SlicedSystem) = A.affine
system(A::SlicedSystem) = A.system

function Base.size(F::SlicedSystem)
    m, n = size(F.system)
    (m + codim(F.affine), n)
end

(F::SlicedSystem)(x::AbstractVector, p = nothing) = F.sliced(x, p)
ModelKit.evaluate!(u, F::SlicedSystem, x::AbstractVector, p = nothing) =
    evaluate!(u, F.sliced, x, p)
ModelKit.evaluate_and_jacobian!(u, U, F::SlicedSystem, x, p = nothing) =
    evaluate_and_jacobian!(u, U, F.sliced, x, p)
ModelKit.taylor!(u, v::Val, F::SlicedSystem, tx, p = nothing) =
    taylor!(u, v, F.sliced, tx, p)