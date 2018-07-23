export OrthogonalPatch

"""
    OrthogonalPatch()

"""
struct OrthogonalPatch <: AbstractLocalAffinePatch end

struct OrthogonalPatchState{T, H} <: AbstractAffinePatchState
    v_conj::PVector{T, H} # this is already conjugated
end

function state(::OrthogonalPatch, x::PVector)
    v = copy(x)
    conj!(raw(v))
    OrthogonalPatchState(v)
end
nequations(::OrthogonalPatchState) = 1

function precondition!(state::OrthogonalPatchState, x)
    LinearAlgebra.normalize!(x)
    update!(state, x)
end
function update!(state::OrthogonalPatchState, x)
    raw(state.v_conj) .= conj.(raw(x))
    LinearAlgebra.normalize!(raw(state.v_conj))
    nothing
end

function evaluate!(u, state::OrthogonalPatchState, x)
    out = -one(eltype(x))
    @inbounds for i=1:length(x)
        out += state.v_conj[i] * x[i]
    end
    u[end] = out
    nothing
end

function jacobian!(U, state::OrthogonalPatchState, x)
    @inbounds for j=1:size(U, 2)
        U[end, j] = state.v_conj[j]
    end
    nothing
end
