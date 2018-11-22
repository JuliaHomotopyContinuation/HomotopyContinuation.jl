import LinearAlgebra

export OrthogonalPatch

"""
    OrthogonalPatch()

"""
struct OrthogonalPatch <: AbstractLocalAffinePatch end

struct OrthogonalPatchState{T, H} <: AbstractAffinePatchState
    v̄::PVector{T, H} # this is already conjugated
end

function state(::OrthogonalPatch, x::PVector)
    v = copy(x)
    conj!(raw(v))
    OrthogonalPatchState(v)
end
nequations(::OrthogonalPatchState) = 1


function precondition!(state::OrthogonalPatchState, x)
    update!(state, x)
end


function normalize!(x::AbstractVector, ::OrthogonalPatchState)
    LinearAlgebra.normalize!(x)
end

Base.@propagate_inbounds function update!(state::OrthogonalPatchState, x, isnormalized=false)
    @boundscheck length(x) == length(state.v̄)
    @inbounds for i in eachindex(state.v̄)
        state.v̄[i] = conj(x[i])
    end
    !isnormalized && LinearAlgebra.normalize!(state.v̄)
    state
end

Base.@propagate_inbounds function evaluate!(u, state::OrthogonalPatchState, x)
    @boundscheck length(x) == length(state.v̄)
    out = -one(eltype(x))
    @inbounds for i in eachindex(state.v̄)
        out += state.v̄[i] * x[i]
    end
    u[end] = out
    u
end

function jacobian!(U, state::OrthogonalPatchState, x)
    @inbounds for j=1:size(U, 2)
        U[end, j] = state.v̄[j]
    end
    nothing
end
