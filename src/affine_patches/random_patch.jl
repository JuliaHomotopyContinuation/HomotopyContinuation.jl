export RandomPatch

"""
    RandomPatch()

A random patch. The vector has norm 1.
"""
struct RandomPatch <: AbstractAffinePatch end

struct RandomPatchState{T, N} <: AbstractAffinePatchState{N}
    v̄::PVector{T, N}
end

is_global_patch(::RandomPatch) = true

function state(::RandomPatch, x::PVector)
    v = copy(x)
    Random.randn!(v)
    LinearAlgebra.normalize!(v)
    conj!(v)
    RandomPatchState(v)
end

onpatch!(x::AbstractVector, state::RandomPatchState) = onpatch!(x, state.v̄)
evaluate!(u, state::RandomPatchState, x::PVector) = evaluate_patch!(u, state.v̄, x)
jacobian!(U, state::RandomPatchState, x::PVector) = jacobian_patch!(U, state.v̄, x)
