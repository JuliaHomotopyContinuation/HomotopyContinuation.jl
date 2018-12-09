export FixedPatch

"""
    FixedPatch()
"""
struct FixedPatch <: AbstractLocalAffinePatch end

struct FixedPatchState{T, H} <: AbstractAffinePatchState
    v̄::PVector{T, H} # this is already conjugated
end

function state(::FixedPatch, x::PVector)
    v = copy(x)
    conj!(raw(v))
    FixedPatchState(v)
end
nequations(::FixedPatchState) = 1

function setup!(state::FixedPatchState, x::AbstractVector)
    @boundscheck length(x) == length(state.v̄)
    LinearAlgebra.normalize!(x)
    @inbounds for i in eachindex(state.v̄)
        state.v̄[i] = conj(x[i])
    end
    state
end

function onpatch!(x::AbstractVector, state::FixedPatchState)
    @boundscheck length(x) == length(state.v̄)
    λ = zero(eltype(x))
    @inbounds for i=1:length(x)
        λ += state.v̄[i] * x[i]
    end
    λ⁻¹ = @fastmath inv(λ)
    for i=1:length(x)
        x[i] *= λ⁻¹
    end
    x
end

function evaluate!(u, state::FixedPatchState, x::PVector)
    out = -one(eltype(x))
    @inbounds for i=1:length(x)
        out += state.v̄[i] * x[i]
    end
    @inbounds u[end] = out
    nothing
end

function jacobian!(U, state::FixedPatchState, x::PVector)
    @inbounds for j=1:size(U, 2)
        U[end, j] = state.v̄[j]
    end
    nothing
end
