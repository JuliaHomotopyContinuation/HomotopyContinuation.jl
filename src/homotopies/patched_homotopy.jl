export PatchedHomotopy, PatchedHomotopyCache

"""
    PatchedHomotopy(H::AbstractHomotopy, patch, v::PVector)

Augment the homotopy `H` with the given patch `v`. This results in the system `[H(x,t); v ⋅ x - 1]`
"""
struct PatchedHomotopy{H<:AbstractHomotopy, PS<:AbstractAffinePatchState} <: AbstractHomotopy
    homotopy::H
    patch::PS
end

function PatchedHomotopy(H::AbstractHomotopy, p::AbstractAffinePatch, x)
    PatchedHomotopy(H, state(p, x))
end

function Base.size(H::PatchedHomotopy)
    m, n = size(H.homotopy)
    (m + nequations(H.patch), n)
end

struct PatchedHomotopyCache{HC, T} <: AbstractHomotopyCache
    cache::HC
    A::Matrix{T} # intermediate storage of the jacobian
    b::Vector{T} # intermediate storage for the evaluation
end

function cache(ph::PatchedHomotopy, x, t)
    H = ph.homotopy
    c = cache(H, x, t)
    PatchedHomotopyCache(c, jacobian(H, x, t, c), evaluate(H, x, t, c))
end

function evaluate!(u, H::PatchedHomotopy, x, t, c::PatchedHomotopyCache)
    M, N = size(H.homotopy)

    # H(x, t)
    evaluate!(c.b, H.homotopy, x, t, c.cache)
    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    evaluate!(u, H.patch, x)
    u
end
function jacobian!(U, H::PatchedHomotopy, x::PVector, t, c::PatchedHomotopyCache)
    M, N = size(H.homotopy)
    # J_H(x, t)
    jacobian!(c.A, H.homotopy, x, t, c.cache)
    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end
    jacobian!(U, H.patch, x)
    U
end

function dt!(u, H::PatchedHomotopy, x::PVector, t, c::PatchedHomotopyCache)
    M, N = size(H.homotopy)
    # [H(x,t); v ⋅ x - 1]/∂t = [∂H(x,t)/∂t; 0]
    dt!(c.b, H.homotopy, x, t, c.cache)
    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    for i=1:nequations(H.patch)
        u[M+i] = zero(eltype(u))
    end
    u
end

function evaluate_and_jacobian!(u, U, H::PatchedHomotopy, x::PVector, t, c::PatchedHomotopyCache)
    M, N = size(H.homotopy)
    evaluate_and_jacobian!(c.b, c.A, H.homotopy, x, t, c.cache)

    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end
    jacobian!(U, H.patch, x)

    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    evaluate!(u, H.patch, x)

    nothing
end

function jacobian_and_dt!(U, u, H::PatchedHomotopy, x::PVector, t, c::PatchedHomotopyCache)
    M, N = size(H.homotopy)
    A, b = c.A, c.b
    jacobian_and_dt!(A, b, H.homotopy, x, t, c.cache)
    # jacobian
    @inbounds for j=1:N, i=1:M
        U[i, j] = A[i, j]
    end
    jacobian!(U, H.patch, x)

    # dt
    @inbounds for i=1:M
        u[i] = b[i]
    end
    for i=1:nequations(H.patch)
        u[M+i] = zero(eltype(u))
    end

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

basehomotopy(H::PatchedHomotopy) = basehomotopy(H.homotopy)
