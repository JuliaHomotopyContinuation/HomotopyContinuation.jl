export AffineSlicedSystem, slice

"""
    AffineSlicedSystem(F::AbstractSystem, L::AffineSubspace)

Given a polynomial system `F` with the affine linear subspace `L`
defined by ``Ax = b`` this constructs the system `[F(x); A x - b] = 0`.
See also [`slice`](@ref).
"""
struct AffineSlicedSystem{S<:AbstractSystem,T} <: AbstractSystem
    system::S
    affine::AffineSubspace{T}
end
AffineSlicedSystem(F::System, A::AffineSubspace; compile::Bool = COMPILE_DEFAULT[]) =
    AffineSlicedSystem(fixed(F; compile = compile), A)

"""
    slice(F::System, L::AffineSubspace; compile = $(COMPILE_DEFAULT[]))
    slice(F::AbstractSystem, L::AffineSubspace)

Slice the zero set of the polynomial system `F` with the affine linear subspace `L`
defined by ``Ax = b``. This constructs the system `[F(x); A x - b] = 0`.
See also [`AffineSlicedSystem`](@ref).
"""
slice(F::System, A::AffineSubspace; compile::Bool = COMPILE_DEFAULT[]) =
    AffineSlicedSystem(fixed(F; compile = compile), A)
slice(F::AbstractSystem, A::AffineSubspace) = AffineSlicedSystem(F, A)
Base.intersect(F::Union{System,AbstractSystem}, A::AffineSubspace) =
    AffineSlicedSystem(F, A)
Base.intersect(A::AffineSubspace, F::Union{System,AbstractSystem}) =
    AffineSlicedSystem(F, A)

ModelKit.variables(F::AffineSlicedSystem) = variables(F.system)
ModelKit.parameters(F::AffineSlicedSystem) = parameters(F.system)
ModelKit.variable_groups(F::AffineSlicedSystem) = nothing

affine_subspace(A::AffineSlicedSystem) = A.affine
system(A::AffineSlicedSystem) = A.system

function Base.size(F::AffineSlicedSystem)
    m, n = size(F.system)
    (m + codim(F.affine), n)
end

function (F::AffineSlicedSystem)(x::AbstractVector, p = nothing)
    u = F.system(x, p)
    append!(u, extrinsic(F.affine)(x))
    u
end

function ModelKit.evaluate!(u, F::AffineSlicedSystem, x, p = nothing)
    evaluate!(u, F.system, x, p)

    @unpack A, b = extrinsic(F.affine)
    m = size(F.system, 1)
    n, N = size(A)
    @inbounds for i = 1:n
        u[m+i] = -b[i]
    end
    @inbounds for j = 1:N
        xj = x[j]
        for i = 1:n
            u[m+i] += A[i, j] * xj
        end
    end

    u
end

function ModelKit.evaluate_and_jacobian!(u, U, F::AffineSlicedSystem, x, p = nothing)
    evaluate_and_jacobian!(u, U, F.system, x, p)

    @unpack A, b = extrinsic(F.affine)
    m = size(F.system, 1)
    n, N = size(A)
    @inbounds for i = 1:n
        u[m+i] = -b[i]
    end
    @inbounds for j = 1:N
        xj = x[j]
        for i = 1:n
            u[m+i] += A[i, j] * xj
        end
    end

    @inbounds for j = 1:N, i = 1:n
        U[m+i, j] = A[i, j]
    end

    nothing
end

function ModelKit.taylor!(u, v::Val, F::AffineSlicedSystem, tx, p = nothing)
    u .= zero(eltype(u))
    taylor!(u, v, F.system, tx, p)
    u
end
