module AffinePatches

import ..ProjectiveVectors: PVector, raw

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

function precondition!(v, x, ::OrthogonalPatch)
    normalize!(x)
    v .= x
end

init_patch(::OrthogonalPatch, x) = x
function update_patch!(v, ::OrthogonalPatch, x)
    v .= x
    nothing
end

end
