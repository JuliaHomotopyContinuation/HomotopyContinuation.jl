export EmbeddingPatch

"""
    EmbeddingPatch()

Holds an `PVector` onto its affine patch. With this the effect
is basically the same as tracking in affine space.
"""
struct EmbeddingPatch <: AbstractAffinePatch end

struct EmbeddingPatchState <: AbstractAffinePatchState
    nequations::Int
end

function state(::EmbeddingPatch, x::PVector{T,N}) where {T,N}
    EmbeddingPatchState(N)
end
nequations(state::EmbeddingPatchState) = state.nequations

function onpatch!(x::PVector, state::EmbeddingPatchState)
    ProjectiveVectors.affine!(x)
end

function evaluate!(u, state::EmbeddingPatchState, x::PVector{T, N}) where {T, N}
    homvars = ProjectiveVectors.homvars(x)
    n = length(u) - N
    for i in 1:N
        u[n + i] = x[homvars[i]] - one(eltype(x))
    end
    nothing
end

function jacobian!(U, state::EmbeddingPatchState, x::PVector{T, N}) where {T, N}
    homvars = ProjectiveVectors.homvars(x)
    n = length(u) - N
    for j in 1:size(U, 2), i in 1:N
        U[n + i, j] = j == homvars[i] ? one(eltype(x)) : zero(eltype(x))
    end
    nothing
end
