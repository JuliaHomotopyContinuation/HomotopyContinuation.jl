export PatchedSystem, PatchedSystemCache

"""
    PatchedSystem(F::AbstractSystem, patch, v::PVector)

Augment the system `F` with the given patch `v`. This results in the system `[F(x,t); v ⋅ x - 1]`
"""
struct PatchedSystem{F<:AbstractSystem, PS<:AbstractAffinePatchState} <: AbstractSystem
    system::F
    patch::PS
end

function PatchedSystem(F::AbstractSystem, p::AbstractAffinePatch, x)
    PatchedSystem(F, state(p, x))
end

function Base.size(F::PatchedSystem)
    m, n = size(F.system)
    (m + nequations(F.patch), n)
end

struct PatchedSystemCache{FC, T} <: AbstractSystemCache
    cache::FC
    A::Matrix{T} # intermediate storage of the jacobian
    b::Vector{T} # intermediate storage for the evaluation
end

function cache(ph::PatchedSystem, x)
    F = ph.system
    c = cache(F, x)
    b = evaluate(F, x, c)
    A = similar(b, size(F))
    PatchedSystemCache(c, A, b)
end

function evaluate!(u, F::PatchedSystem, x, c::PatchedSystemCache)
    M, N = size(F.system)

    # F(x)
    evaluate!(c.b, F.system, x, c.cache)
    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    evaluate!(u, F.patch, x)
    u
end

function jacobian!(U, F::PatchedSystem, x::PVector, c::PatchedSystemCache)
    M, N = size(F.system)
    # J_F(x)
    jacobian!(c.A, F.system, x, c.cache)
    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end
    jacobian!(U, F.patch, x)
    U
end

function evaluate_and_jacobian!(u, U, F::PatchedSystem, x::PVector, c::PatchedSystemCache)
    M, N = size(F.system)
    evaluate_and_jacobian!(c.b, c.A, F.system, x, c.cache)

    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end
    jacobian!(U, F.patch, x)

    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    evaluate!(u, F.patch, x)

    nothing
end

function evaluate(F::PatchedSystem, x, c::PatchedSystemCache)
    M, N = size(F)
    evaluate!(similar(c.b, M), F, x, c)
end

function jacobian(F::PatchedSystem, x, c::PatchedSystemCache)
    M, N = size(F)
    jacobian!(similar(c.A, M, N), F, x, c)
end
