module AffinePatches

import ..ProjectiveVectors
import ..ProjectiveVectors: AbstractProjectiveVector, PVector, raw

export AbstractAffinePatch,
    AbstractLocalAffinePatch,
    state,
    precondition!,
    update!

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
    state(::AbstractAffinePatch, x)::AbstractAffinePatchState

Construct the state of the path from `x`.
"""
function state end

"""
    precondition!(v::AbstractAffinePatchState, x)

Modify both such that `v` is properly setup and `v⋅x-1=0` holds.
"""
function precondition! end

"""
    update_patch!(::AbstractAffinePatchState, x)

Update the patch depending on the local state.
"""
update!(::AbstractAffinePatchState, x) = nothing


include("affine_patches/orthogonal_patch.jl")
include("affine_patches/embedding_patch.jl")
include("affine_patches/random_patch.jl")
include("affine_patches/fixed_patch.jl")

end
