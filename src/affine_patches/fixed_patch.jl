export FixedPatch

"""
    FixedPatch()

"""
struct FixedPatch <: AbstractLocalAffinePatch end

struct FixedPatchState{V<:AbstractProjectiveVector} <: AbstractAffinePatchState
    v_conj::V # this is already conjugated
end

function state(::FixedPatch, x::AbstractProjectiveVector)
    v = copy(x)
    conj!(raw(v))
    FixedPatchState(v)
end
nequations(::FixedPatchState{<:PVector}) = 1

function precondition!(state::FixedPatchState, x)
    raw(state.v_conj) .= conj.(raw(x)) ./ sum(abs2, x)
    nothing
end

function evaluate!(u, state::FixedPatchState{<:PVector}, x)
    out = -one(eltype(x))
    @inbounds for i=1:length(x)
        out += state.v_conj[i] * x[i]
    end
    u[end] = out
    nothing
end

function jacobian!(U, state::FixedPatchState{<:PVector}, x)
    for j=1:size(U, 2)
        U[end, j] = state.v_conj[j]
    end
    nothing
end
