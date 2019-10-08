export AbstractAffinePatch, nequations, state, onpatch!, setup!, changepatch!

"""
    AbstractAffinePatch

An affine patch is a hyperplane defined by ``v⋅x-1=0``. The currently implemented affine patches
are:
* [`EmbeddingPatch`](@ref)
* [`FixedPatch`](@ref)
* [`OrthogonalPatch`](@ref)
* [`RandomPatch`](@ref)
"""
abstract type AbstractAffinePatch end

"""
    AbstractAffinePatchState{N}

This holds the actual patch information. `N` indicates the number of products of projective spaces.
"""
abstract type AbstractAffinePatchState{N} end

"""
    nequations(::AbstractAffinePatchState)

Number of equations an affine patch adds.
"""
nequations(::AbstractAffinePatchState{N}) where {N} = N

"""
    state(::AbstractAffinePatch, x)::AbstractAffinePatchState

Construct the state of the path from `x`.
"""
state(p::AbstractAffinePatch, x) = throw(MethodError(state, (p, x)))

"""
    onpatch!(x::AbstractVector, ::AbstractAffinePatchState)

Scale a vector `x` such that it is on the affine patch.
"""
function onpatch! end

"""
    init!(::AbstractAffinePatchState, x::AbstractVector)

Setup the affine patch depending on `x` and modify `x` if necessary.
This is only called once at the beginning of a tracked path.
"""
init!(state::AbstractAffinePatchState, x::AbstractVector) = onpatch!(x, state)
@deprecate setup!(state::AbstractAffinePatchState, x::AbstractVector) init!(state, x)


"""
    changepatch!(::AbstractAffinePatch, x::AbstractVector)

The same as [`setup!`](@ref) but only called during the path tracking.
"""
changepatch!(::AbstractAffinePatchState, x::AbstractVector) = nothing

is_global_patch(::AbstractAffinePatch) = false

include("affine_patches/common.jl")
include("affine_patches/orthogonal_patch.jl")
include("affine_patches/embedding_patch.jl")
include("affine_patches/random_patch.jl")
include("affine_patches/fixed_patch.jl")
