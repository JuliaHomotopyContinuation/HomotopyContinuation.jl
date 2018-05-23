import Base: isless, collect
import ..ProjectiveVectors
import ..Utilities

## The structure for sorting the solutions is a binary search tree (BST).
## First, we cluster the points by ordering them w.r.t the the absolute value of their first entry.
## Then, we compare the points pairwise in each cluster.

mutable struct BST
    pos::Vector{Int}
    left::Nullable{BST}
    right::Nullable{BST}
end

function Base.isless(x::Vector{T}, y::Vector{S}) where {S,T <: Number}
    if abs(normalize(x)[1]) < abs(normalize(y)[1])
        return true
    else
        return false
    end
end

function compare(x, y, tol)
    # if ProjectiveVectors.infinity_norm(x, y) > tol
    if infinity_norm(x-y) > tol
        return false
    else
        return true
    end
end



function push_for_clustering!(node::BST, i, vectors, τ)
    if compare(vectors[i], vectors[node.pos[1]], τ)
        push!(node.pos, i)
    else
        if vectors[i] < vectors[node.pos[1]]
            if isnull(node.left)
                node.left = BST([i], Nullable{BST}(), Nullable{BST}())
            else
                push_for_clustering!(node.left.value, i, vectors, τ)
            end
        else
            if isnull(node.right)
                node.right = BST([i], Nullable{BST}(), Nullable{BST}())
            else
                push_for_clustering!(node.right.value, i, vectors, τ)
            end
        end
    end
end

function push_for_identifying_multiplicities!(node::BST, i, vectors, τ)
    if compare(vectors[i], vectors[node.pos[1]], τ)
        push!(node.pos, i)
    else
        if isnull(node.left)
            node.left = BST([i], Nullable{BST}(), Nullable{BST}())
        else
            push_for_identifying_multiplicities!(node.left.value, i, vectors, τ)
        end
    end
end

function Base.collect(node::BST)
    if isnull(node.left)
        if isnull(node.right)
            return [node.pos]
        else
            return [[node.pos]; collect(node.right.value)]
        end
    else
        if isnull(node.right)
            return [[node.pos]; collect(node.left.value)]
        else
            return [[node.pos]; collect(node.left.value); collect(node.right.value)]
        end
    end
end


"""
    multiplicities(V, tol)

Returns an array of arrays of integers. Each array v in V contains all indices i,j such that V[i] and V[j] have distance at most tol.
"""

function multiplicities(V::Vector{Vector{T}}, tol) where {T<:Number}

    root = BST([1], Nullable{BST}(), Nullable{BST}())
    for i in 2:length(V)
        push_for_clustering!(root, i, V, tol)
    end

    M = collect(root)
    mults = Vector{Int}[]

    for m in M
        root = BST([1], Nullable{BST}(), Nullable{BST}())
        for i in 2:length(m)
            push_for_identifying_multiplicities!(root, i, V[m], tol)
        end
        for n in collect(root)
            push!(mults, m[n])
        end
    end

    mults
end
