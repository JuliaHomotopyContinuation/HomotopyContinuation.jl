export FixedPatch

"""
    FixedPatch()
"""
struct FixedPatch <: AbstractLocalAffinePatch end

struct FixedPatchState{T, N} <: AbstractAffinePatchState
    v̄::PVector{T, N} # this is already conjugated
end

function state(::FixedPatch, x::PVector)
    v = copy(x)
    conj!(v.data)
    FixedPatchState(v)
end
nequations(::FixedPatchState{T, N}) where {T, N}= N

function setup!(state::FixedPatchState, x::AbstractVector)
    @boundscheck length(x) == length(state.v̄)
    LinearAlgebra.normalize!(x)
    @inbounds for i in eachindex(state.v̄)
        state.v̄[i] = conj(x[i])
    end
    state
end

onpatch!(x::AbstractVector, state::FixedPatchState) = onpatch!(x, state.v̄)
evaluate!(u, state::FixedPatchState, x::PVector) = evaluate!(u, state.v̄, x)
jacobian!(U, state::FixedPatchState, x::PVector) = jacobian!(U, state.v̄, x)
