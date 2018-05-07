export PatchSwitcherHomotopy, PatchSwitcherHomotopyCache

import ..AffinePatches
import ..AffinePatches: AbstractAffinePatch, AbstractAffinePatchState
import ..ProjectiveVectors: AbstractProjectiveVector, PVector, raw

"""
    PatchSwitcherHomotopy(H::AbstractHomotopy, patch, v::AbstractProjectiveVector)

Augment the homotopy `H` with the given patch `v`. This results in the system `[H(x,t); v ⋅ x - 1]`
"""
struct PatchSwitcherHomotopy{H<:Systems.FixedHomotopy, P1<:AbstractAffinePatchState, P2<:AbstractAffinePatchState} <: AbstractHomotopy
    homotopy::H
    start_patch::P1
    target_patch::P2
end

function Base.size(H::PatchSwitcherHomotopy)
    m, n = size(H.homotopy)
    (m + AffinePatches.nequations(H.start_patch), n)
end

struct PatchSwitcherHomotopyCache{HC, T} <: AbstractHomotopyCache
    cache::HC
    A::Matrix{T} # intermediate storage of the jacobian
    b::Vector{T} # intermediate storage for the evaluation
    patch_jac::Matrix{T}
    patch_val::Vector{T}
end

function cache(ph::PatchSwitcherHomotopy, x, t)
    H = ph.homotopy
    _, n = size(H)
    m = AffinePatches.nequations(ph.start_patch)

    c = Systems.cache(H, x)
    A, b = Systems.jacobian(H, x, c), Systems.evaluate(H, x, c)
    patch_jac = similar(A, m, n)
    patch_val = similar(b, m)
    PatchSwitcherHomotopyCache(c, A, b, patch_jac, patch_val)
end

precondition!(H::PatchSwitcherHomotopy, x, t, cache) = AffinePatches.precondition!(H.start_patch, x)

function evaluate!(u, H::PatchSwitcherHomotopy, x, t, c::PatchSwitcherHomotopyCache)
    M, _ = size(H.homotopy)
    # H(x, t)
    Systems.evaluate!(c.b, H.homotopy, raw(x), c.cache)
    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    _patches_evaluate!(u, H, x, t, c)
    u
end

function _patches_evaluate!(u, H, x, t, c)
    M, N = size(H.homotopy)
    m = AffinePatches.nequations(H.start_patch)
    AffinePatches.evaluate!(c.patch_val, H.start_patch, x)
    AffinePatches.evaluate!(u, H.target_patch, x)
    for i=1:m
        u[M+i] = (1-t) * u[M+i] + t * c.patch_val[i]
    end
    nothing
end

function jacobian!(U, H::PatchSwitcherHomotopy, x::AbstractProjectiveVector, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H.homotopy)
    # J_H(x, t)
    Systems.jacobian!(c.A, H.homotopy, raw(x), c.cache)
    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end

    _patches_jacobian!(U, H, x, t, c)

    U
end

function _patches_jacobian!(U, H, x, t, c)
    M, N = size(H.homotopy)
    m = AffinePatches.nequations(H.start_patch)
    AffinePatches.jacobian!(c.patch_jac, H.start_patch, x)
    AffinePatches.jacobian!(U, H.target_patch, x)
    for j=1:N, i=1:m
        U[M+i, j] = (1-t) * U[M+i, j] + t * c.patch_jac[i, j]
    end
    nothing
end

function dt!(u, H::PatchSwitcherHomotopy, x::AbstractProjectiveVector, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H.homotopy)
    # [H(x,t); v ⋅ x - 1]/∂t = [∂H(x,t)/∂t; 0]
    @inbounds for i=1:M
        u[i] = zero(eltype(u))
    end
    for i=1:AffinePatches.nequations(H.start_patch)
        u[M+i] = zero(eltype(u))
    end
    u
end

function _patches_dt!(u, H, x, t, c)
    M, _ = size(H.homotopy)
    m = AffinePatches.nequations(H.start_patch)
    AffinePatches.evaluate!(c.patch_val, H.start_patch, x)
    AffinePatches.evaluate!(u, H.target_patch, x)
    for i=1:m
        u[M+i] = c.patch_val[i] - u[M+i]
    end
    nothing
end


function evaluate_and_jacobian!(u, U, H::PatchSwitcherHomotopy, x::AbstractProjectiveVector, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H.homotopy)
    Systems.evaluate_and_jacobian!(c.b, c.A, H.homotopy, raw(x), c.cache)

    @inbounds for j=1:N, i=1:M
        U[i, j] = c.A[i, j]
    end
    _patches_jacobian!(U, H, x, t, c)

    @inbounds for i=1:M
        u[i] = c.b[i]
    end
    _patches_evaluate!(U, H, x, t, c)

    nothing
end

function evaluate(H::PatchSwitcherHomotopy, x, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H)
    evaluate!(similar(c.b, M), H, x, t, c)
end
function jacobian(H::PatchSwitcherHomotopy, x, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H)
    jacobian!(similar(c.A, M, N), H, x, t, c)
end
function dt(H::PatchSwitcherHomotopy, x, t, c::PatchSwitcherHomotopyCache)
    M, N = size(H)
    dt!(similar(c.b, M), H, x, t, c)
end
