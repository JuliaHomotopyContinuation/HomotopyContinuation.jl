module AffinePatches

import ..ProjectiveVectors: AbstractProjectiveVector, PVector, raw

export AbstractPatch,
    AbstractLocalAffinePatch,
    precondition!,
    init_patch,
    update_patch!,
    OrthogonalPatch

# TODO: This currently does not support products of project spaces. Stil need to figure
# out how to do this best. So I'm punting on that for now...
"""
    AbstractAffinePatch

An affine patch is a hyperplane defined by ``v⋅x-1=0``.
"""
abstract type AbstractAffinePatch end


abstract type AbstractAffinePatchState end

"""
    AbstractLocalAffinePatch

An affine patch is a hyperplane defined by ``v⋅x-1=0``.
"""
abstract type AbstractLocalAffinePatch <: AbstractAffinePatch end


"""
    precondition!(v, x, ::AbstractAffinePatch)

This is called at the beginning of a tracked path.
`v` is the patch and `x` is the start solution. Modify both such that
`v` is properly setup and `v⋅x-1=0` holds.
"""
function precondition! end


"""
    init_patch(::AbstractAffinePatch, x)

Initialize the patch.
"""
function init_patch end


"""
    update_patch!(v, ::AbstractLocalAffinePatch, x)

Update the patch depending on the local state. This is called after
each successfull step.
"""
function update_patch! end

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




#
#
# function update_patch!(v, ::OrthogonalPatch, x)
#     normalize!(x)
#     v .= x
#     nothing
# end

end
