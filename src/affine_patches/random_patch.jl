export RandomPatch


"""
    RandomPatch()

A random patch. The vector has norm 1.
"""
struct RandomPatch <: AbstractAffinePatch end

struct RandomPatchState{T, H} <: AbstractAffinePatchState
    v::PVector{T, H}
end

function state(::RandomPatch, x::AbstractProjectiveVector)
    v = copy(x)
    for i=1:length(v)
        v[i] = rand(eltype(v))
    end
    normalize!(v)
    RandomPatchState(v)
end
nequations(::RandomPatchState) = 1

function precondition!(state::RandomPatchState, x::PVector)
    scale!(raw(x), inv(dot(state.v, x)))
end

function evaluate!(u, state::RandomPatchState, x)
    u[end] = dot(state.v, x) - one(eltype(x))
    nothing
end

function jacobian!(U, state::RandomPatchState, x)
    for j=1:size(U, 2)
        U[end, j] = conj(state.v[j])
    end
    nothing
end
