# import ..AffinePatches: AbstractAffinePatch
import ..Homotopies: AbstractHomotopy, AbstractHomotopyCache, HomotopyWithCache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt
import ..Homotopies
import ..ProjectiveVectors: AbstractProjectiveVector, PVector, raw

"""
    PatchedHomotopy(H::AbstractHomotopy, v::AbstractProjectiveVector)

Augment the homotopy `H` with the given patch `v`. This results in the system `[H(x,t); v ⋅ x - 1]`
"""
struct PatchedHomotopy{M, N, H<:AbstractHomotopy, V<:AbstractProjectiveVector} <: AbstractHomotopy{M, N}
    homotopy::H
    patch::V
end

function PatchedHomotopy(hom::H, patch::V) where {M, N, H<:AbstractHomotopy{M, N}, V<:PVector}
   PatchedHomotopy{M+1, N, H, V}(hom, patch)
end

"""
    patch(H::PatchedHomotopy)

Get the used patch.
"""
patch(H::PatchedHomotopy) = H.patch

struct PatchedHomotopyCache{HC, T} <: AbstractHomotopyCache
    cache::HC
    A::Matrix{T} # intermediate storage of the jacobian
    b::Vector{T} # intermediate storage for the evaluation
end

function Homotopies.cache(ph::PatchedHomotopy, x, t)
    H = HomotopyWithCache(ph.homotopy, x, t)
    PatchedHomotopyCache(H.cache, jacobian(H, x, t), H(x, t))
end

function evaluate!(u, H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    # H(x, t)
    evaluate!(c.b, H.homotopy, x, t, c.cache)
    @inbounds for i=1:(M-1)
        u[i] = c.b[i]
    end
    # v⋅x - 1
    out = -one(eltype(x))
    @inbounds for i=1:N
        out = muladd(conj(H.patch[i]), x[i], out)
    end
    @inbounds u[M] = out
    u
end
function jacobian!(U, H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    # J_H(x, t)
    jacobian!(c.A, H.homotopy, x, t, c.cache)
    @inbounds for j=1:N, i=1:(M-1)
        U[i, j] = c.A[i, j]
    end
    # gradient of v⋅x - 1 => v'
    @inbounds for j=1:N
        U[M, j] = conj(H.patch[j])
    end
    U
end

function dt!(u, H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    # [H(x,t); v ⋅ x - 1]/∂t = [∂H(x,t)/∂t; 0]
    dt!(c.b, H.homotopy, x, t, c.cache)
    @inbounds for i=1:(M-1)
        u[i] = c.b[i]
    end
    @inbounds u[M] = zero(eltype(u))
    u
end

function evaluate_and_jacobian!(u, U, H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    evaluate_and_jacobian!(c.b, c.A, H.homotopy, x, t, c.cache)
    # jacobian
    @inbounds for j=1:N, i=1:(M-1)
        U[i, j] = c.A[i, j]
    end
    @inbounds for j=1:N
        U[M, j] = conj(H.patch[j])
    end
    # eval
    @inbounds for i=1:(M-1)
        u[i] = c.b[i]
    end
    # v⋅x - 1
    out = -one(eltype(x))
    @inbounds for i=1:N
        out = muladd(conj(H.patch[i]), x[i], out)
    end
    @inbounds u[M] = out

    nothing
end

function jacobian_and_dt!(U, u, H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    A, b = c.A, c.b
    jacobian_and_dt!(A, b, H.homotopy, x, t, c.cache)
    # jacobian
    @inbounds for j=1:N, i=1:(M-1)
        U[i, j] = A[i, j]
    end
    # gradient of v⋅x - 1 => v'
    @inbounds for j=1:N
        U[M, j] = conj(H.patch[j])
    end
    # dt
    @inbounds for i=1:(M-1)
        u[i] = b[i]
    end
    u[M] = zero(eltype(u))

    nothing
end

function evaluate(H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    evaluate!(similar(c.b, M), H, x, t, c)
end
function jacobian(H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    jacobian!(similar(c.A, M, N), H, x, t, c)
end
function dt(H::PatchedHomotopy{M, N}, x, t, c::PatchedHomotopyCache) where {M, N}
    dt!(similar(c.b, M), H, x, t, c)
end
function jacobian_and_dt(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    jacobian_and_dt!(similar(c.A, M, N), similar(c.b, M), H, x, t, c)
end
function evaluate_and_jacobian(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    evaluate_and_jacobian!(similar(c.b, M), similar(c.A, M, N), H, x, t, c)
end
