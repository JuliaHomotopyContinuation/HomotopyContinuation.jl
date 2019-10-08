export OrthogonalPatch

"""
    OrthogonalPatch()

The orthogonal patch is a dynamically changing patch. It computes in such a way that
`||x||=1` with respect to the 2-norm. See (for example) Section 3.2 in [^HR18].

[^HR18]: https://arxiv.org/pdf/1710.06362.pdf
"""
struct OrthogonalPatch <: AbstractAffinePatch end

struct OrthogonalPatchState{T,N} <: AbstractAffinePatchState{N}
    v̄::PVector{T,N} # this is already conjugated
end

function state(::OrthogonalPatch, x::PVector)
    v = copy(x)
    conj!(v.data)
    OrthogonalPatchState(v)
end

function init!(state::OrthogonalPatchState, x::AbstractVector)
    @boundscheck length(x) == length(state.v̄)
    LinearAlgebra.normalize!(x)
    @inbounds for i in eachindex(state.v̄)
        state.v̄[i] = conj(x[i])
    end
    state
end
@propagate_inbounds changepatch!(state::OrthogonalPatchState, x::AbstractVector) =
    init!(state, x)

onpatch!(x::AbstractVector, state::OrthogonalPatchState) = onpatch!(x, state.v̄)
evaluate!(u, state::OrthogonalPatchState, x::PVector) = evaluate_patch!(u, state.v̄, x)
jacobian!(U, state::OrthogonalPatchState, x::PVector) = jacobian_patch!(U, state.v̄, x)
