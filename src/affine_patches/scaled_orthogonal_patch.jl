import LinearAlgebra, ProjectiveVectors
import ..Utilities

export ScaledOrthogonalPatch

"""
    ScaledOrthogonalPatch()

"""
struct ScaledOrthogonalPatch <: AbstractLocalAffinePatch end

mutable struct ScaledOrthogonalPatchState{T, N} <: AbstractAffinePatchState
    v::PVector{Complex{T}, N}
    norm2::NTuple{N, Float64}
end

function state(::ScaledOrthogonalPatch, x::PVector)
    v = copy(x)
    ScaledOrthogonalPatchState(v, LinearAlgebra.norm(v))
end
nequations(::ScaledOrthogonalPatchState{T, N}) where {T,N} = N

function setup!(state::ScaledOrthogonalPatchState, x::PVector, ip::InnerProduct)
    @boundscheck length(x) == length(state.v)
    for i in 1:length(x)
        state.v[i] = x[i] / ip[i]^2
    end
    state.norm2 = Utilities.norm2(x, ip)
    state
end
Base.@propagate_inbounds function changepatch!(state::ScaledOrthogonalPatchState, x::PVector, ip::InnerProduct)
    setup!(state, x, ip)
end

function evaluate!(u, state::ScaledOrthogonalPatchState, x::PVector{T,N}) where {T,N}
    v = state.v
    ranges = ProjectiveVectors.dimension_indices(v)
    n = length(u) - N
    for (k, range) in enumerate(ranges)
        out = zero(eltype(x))
        for i in range
            out += conj(v[i]) * x[i]
        end
        u[n + k] = out - state.norm2[k]
    end
    nothing
end

function jacobian!(U, state::ScaledOrthogonalPatchState, x::PVector{T,N}) where {T,N}
    v = state.v
    ranges = ProjectiveVectors.dimension_indices(v)
    n = size(U, 1) - N
    for (k, range) in enumerate(ranges)
        for j in range
            U[n + k, j] = conj(v[j])
        end
    end
    nothing
end
