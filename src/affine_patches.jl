module AffinePatches

import ..ProjectiveVectors
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
abstract type AbstractLocalAffinePatch <: AbstractAffinePatch end

"""
    AbstractAffinePatchState

This holds the actual patch information.
"""
abstract type AbstractAffinePatchState end


"""
    precondition!(v::AbstractAffinePatchState, x)

Modify both such that `v` is properly setup and `v⋅x-1=0` holds.
"""
precondition!


"""
    update_patch!(::AbstractAffinePatchState, x)

Update the patch depending on the local state.
"""
update!(::AbstractAffinePatchState, x) = nothing




include("affine_patches/orthogonal_patch.jl")
include("affine_patches/embedding_patch.jl")
include("affine_patches/random_patch.jl")

end
