export EmbeddingPatch

"""
    EmbeddingPatch()

Holds an `AbstractProjectiveVector` onto its affine patch. With this the effect
is basically the same as tracking in affine space.
"""
struct EmbeddingPatch <: AbstractAffinePatch end

struct EmbeddingPatchState <: AbstractAffinePatchState
    nequations::Int
end

function state(::EmbeddingPatch, x::PVector)
    EmbeddingPatchState(1)
end
nequations(state::EmbeddingPatchState) = state.nequations

function precondition!(state::EmbeddingPatchState, x::AbstractProjectiveVector)
    ProjectiveVectors.affine!(x)
end

function evaluate!(u, state::EmbeddingPatchState, x::PVector{T, Int}) where T
    u[end] = x[ProjectiveVectors.homvar(x)] - one(T)
    nothing
end

function jacobian!(U, state::EmbeddingPatchState, x::PVector{T, Int}) where T
    i = ProjectiveVectors.homvar(x)
    for j=1:size(U, 2)
        U[end, j] = j == i ? one(T) : zero(T)
    end
    nothing
end
