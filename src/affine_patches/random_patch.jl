export RandomPatch

"""
    RandomPatch()

A random patch. For this we first draw entries of a vector `v` independently from a
complex normal distribution (`randn(ComplexF64)`). And then normalize `v` with respect
to the 2-norm.
"""
struct RandomPatch <: AbstractAffinePatch end

struct RandomPatchState{T,N} <: AbstractAffinePatchState{N}
    v::PVector{T,N}
end

is_global_patch(::RandomPatch) = true

function state(::RandomPatch, x::PVector)
    v = similar(x, ComplexF64)
    Random.randn!(v)
    LinearAlgebra.normalize!(v)
    RandomPatchState(v)
end

onpatch!(x::AbstractVector, state::RandomPatchState) = onpatch!(x, state.v)
evaluate!(u, state::RandomPatchState, x::PVector) = evaluate_patch!(u, state.v, x)
jacobian!(U, state::RandomPatchState, x::PVector) = jacobian_patch!(U, state.v, x)
