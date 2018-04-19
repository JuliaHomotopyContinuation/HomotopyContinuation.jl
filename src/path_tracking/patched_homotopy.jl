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
struct PatchedHomotopy{H<:AbstractHomotopy, V<:AbstractProjectiveVector} <: AbstractHomotopy
    homotopy::H
    patch::V
    size::NTuple{2, Int}
end

function PatchedHomotopy(H::AbstractHomotopy, patch::PVector)
    m, n = size(H)
    PatchedHomotopy(H, patch, (m+1, n))
end

Base.size(H::PatchedHomotopy{Hom, <:PVector}) where Hom = H.size

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

function evaluate!(u, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
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
function jacobian!(U, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
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

function dt!(u, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    # [H(x,t); v ⋅ x - 1]/∂t = [∂H(x,t)/∂t; 0]
    dt!(c.b, H.homotopy, x, t, c.cache)
    @inbounds for i=1:(M-1)
        u[i] = c.b[i]
    end
    @inbounds u[M] = zero(eltype(u))
    u
end

function evaluate_and_jacobian!(u, U, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
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

function jacobian_and_dt!(U, u, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
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

function evaluate(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    evaluate!(similar(c.b, M), H, x, t, c)
end
function jacobian(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    jacobian!(similar(c.A, M, N), H, x, t, c)
end
function dt(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    dt!(similar(c.b, M), H, x, t, c)
end
function jacobian_and_dt(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    jacobian_and_dt!(similar(c.A, M, N), similar(c.b, M), H, x, t, c)
end
function evaluate_and_jacobian(H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H)
    evaluate_and_jacobian!(similar(c.b, M), similar(c.A, M, N), H, x, t, c)
end
