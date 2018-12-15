import Random

export RandomPatch


"""
    RandomPatch()

A random patch. The vector has norm 1.
"""
struct RandomPatch <: AbstractAffinePatch end

struct RandomPatchState{T, N} <: AbstractAffinePatchState
    v̄::PVector{T, N}
end

function state(::RandomPatch, x::PVector)
    v = copy(x)
    Random.randn!(v)
    LinearAlgebra.normalize!(v)
    conj!(v)
    RandomPatchState(v)
end
nequations(::RandomPatchState{T,N}) where {T,N} = N

onpatch!(x::AbstractVector, state::RandomPatchState) = onpatch!(x, state.v̄)
evaluate!(u, state::RandomPatchState, x::PVector) = evaluate!(u, state.v̄, x)
jacobian!(U, state::RandomPatchState, x::PVector) = jacobian!(U, state.v̄, x)
