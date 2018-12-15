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

function onpatch!(x::AbstractVector, state::OrthogonalPatchState)
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

function setup!(state::OrthogonalPatchState, x::AbstractVector)
    @boundscheck length(x) == length(state.v̄)
    LinearAlgebra.normalize!(x)
    @inbounds for i in eachindex(state.v̄)
        state.v̄[i] = conj(x[i])
    end
    state
end

Base.@propagate_inbounds changepatch!(state::OrthogonalPatchState, x::AbstractVector) = setup!(state, x)

function evaluate!(u, state::OrthogonalPatchState, x)
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
