module AffinePatches

import ..ProjectiveVectors
import ..ProjectiveVectors: AbstractProjectiveVector, PVector, raw
import LinearAlgebra

export AbstractAffinePatch,
    AbstractLocalAffinePatch,
    state,
    onpatch!,
    setup!,
    changepatch!

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
    onpatch!(x::AbstractVector, ::AbstractAffinePatchState)

Scale a vector `x` such that it is on the affine patch.
"""
function onpatch! end

"""
    setup!(::AbstractAffinePatchState, x::AbstractVector)

Setup the affine patch depending on `x` and modify `x` if necessary.
This is only called once at the beginning of a tracked path.
"""
setup!(state::AbstractAffinePatchState, x::AbstractVector) = onpatch!(x, state)

"""
    changepatch!(::AbstractAffinePatch, x::AbstractVector)

The same as [`setup!`](@ref) but only called during the path tracking.
"""
changepatch!(::AbstractAffinePatchState, x::AbstractVector) = nothing

include("affine_patches/orthogonal_patch.jl")
include("affine_patches/embedding_patch.jl")
include("affine_patches/random_patch.jl")
include("affine_patches/fixed_patch.jl")

end
