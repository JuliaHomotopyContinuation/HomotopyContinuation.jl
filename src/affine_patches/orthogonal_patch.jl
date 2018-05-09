export OrthogonalPatch

"""
    OrthogonalPatch()

"""
struct OrthogonalPatch <: AbstractLocalAffinePatch end

struct OrthogonalPatchState{V<:AbstractProjectiveVector} <: AbstractAffinePatchState
    v_conj::V # this is already conjugated
end

function state(::OrthogonalPatch, x::AbstractProjectiveVector)
    v = copy(x)
    conj!(raw(v))
    OrthogonalPatchState(v)
end
nequations(::OrthogonalPatchState{<:PVector}) = 1

function precondition!(state::OrthogonalPatchState, x)
    normalize!(x)
    update!(state, x)
end
function update!(state::OrthogonalPatchState, x)
    raw(state.v_conj) .= conj.(raw(x))
    nothing
end

function evaluate!(u, state::OrthogonalPatchState{<:PVector}, x)
    out = -one(eltype(x))
    @inbounds for i=1:length(x)
        out += state.v_conj[i] * x[i]
    end
    u[end] = out
    nothing
end

function jacobian!(U, state::OrthogonalPatchState{<:PVector}, x)
    for j=1:size(U, 2)
        U[end, j] = state.v_conj[j]
    end
    nothing
end
