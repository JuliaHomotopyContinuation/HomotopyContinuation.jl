module AffinePatches

import ..NewHomotopies

export AbstractPatch,
    AbstractLocalAffinePatch,
    precondition!,
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
    precondition!(x, ::AbstractAffinePatch)

Precondition the start solution `x` for the use of this affine patch.
"""
function precondition! end

"""
    update_patch!(v, ::AbstractLocalAffinePatch, H, x, t)

Update the patch depending on the local state. This is called after
each successfull step.
"""
function update_patch! end

struct OrthogonalPatch <: AbstractLocalAffinePatch end

precondition!(x, ::OrthogonalPatch) = normalize!(x)
function update_patch!(v, ::OrthogonalPatch, H, x, t)
    v .= x
    nothing
end

end
