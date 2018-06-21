export FixedPatch

"""
    FixedPatch()
"""
struct FixedPatch <: AbstractLocalAffinePatch end

struct FixedPatchState{T, H} <: AbstractAffinePatchState
    v_conj::PVector{T, H} # this is already conjugated
end

function state(::FixedPatch, x::PVector)
    v = copy(x)
    conj!(raw(v))
    FixedPatchState(v)
end
nequations(::FixedPatchState) = 1

function precondition!(state::FixedPatchState, x::PVector)
    raw(state.v_conj) .= conj.(raw(x)) ./ sum(abs2, x)
    nothing
end

function evaluate!(u, state::FixedPatchState, x::PVector)
    out = -one(eltype(x))
    @inbounds for i=1:length(x)
        out += state.v_conj[i] * x[i]
    end
    u[end] = out
    nothing
end

function jacobian!(U, state::FixedPatchState, x::PVector)
    for j=1:size(U, 2)
        U[end, j] = state.v_conj[j]
    end
    nothing
end
