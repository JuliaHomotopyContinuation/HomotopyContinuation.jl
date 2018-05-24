using ..Utilities

## The structure for sorting the solutions is a binary search tree (BST).
## First, we cluster the points by ordering them w.r.t the the absolute value of their first entry.
## Then, we compare the points pairwise in each cluster.
mutable struct BST
    pos::Vector{Int}
    left::Nullable{BST}
    right::Nullable{BST}
end

function push_for_clustering!(node::BST, i, vectors, τ)
    xᵢ, y = abs(vectors[i][1]), abs(vectors[node.pos[1]][1])
    if abs(xᵢ - y) < τ
        push!(node.pos, i)
    elseif xᵢ < y
        if isnull(node.left)
            node.left = BST([i], Nullable{BST}(), Nullable{BST}())
        else
            push_for_clustering!(get(node.left), i, vectors, τ)
        end
    else
        if isnull(node.right)
            node.right = BST([i], Nullable{BST}(), Nullable{BST}())
        else
            push_for_clustering!(get(node.right), i, vectors, τ)
        end
    end
end

function push_for_identifying_multiplicities!(node::BST, i, vectors, τ)
    # This compares with the infinity_norm but avoids a temporary vector allocation
    if infinity_norm(vectors[i], vectors[node.pos[1]]) < τ
        push!(node.pos, i)
    else
        if isnull(node.left)
            node.left = BST([i], Nullable{BST}(), Nullable{BST}())
        else
            push_for_identifying_multiplicities!(node.left.value, i, vectors, τ)
        end
    end
end

function Base.foreach(f::F, node::BST) where {F<:Function}
    f(node.pos)
    if !isnull(node.left)
        foreach(f, get(node.left))
    end
    if !isnull(node.right)
        foreach(f, get(node.right))
    end
    nothing
end

# Currently we do not need this since `foreach` is sufficient, but maybe its useful later
#
# Base.collect(node::BST) = collectpositions!(Vector{Int}[], node)
# function collectpositions!(positions, node::BST)
#     push!(positions, node.pos)
#     if !isnull(node.left)
#         collectpositions!(positions, get(node.left))
#     end
#     if !isnull(node.right)
#         collectpositions!(positions, get(node.right))
#     end
#     positions
# end
#

"""
    multiplicities(V, tol)

Returns an array of arrays of integers. Each array v in V contains all indices i,j such that V[i] and V[j] have distance at most tol.
"""
function multiplicities(V::Vector{Vector{T}}, tol) where {T<:Number}
    root = BST([1], Nullable{BST}(), Nullable{BST}())
    for i in 2:length(V)
        push_for_clustering!(root, i, V, tol)
    end

    mults = Vector{Int}[]
    foreach(root) do m
        clusterroot = BST([1], Nullable{BST}(), Nullable{BST}())
        for i in 2:length(m)
            push_for_identifying_multiplicities!(clusterroot, i, V[m], tol)
        end
        foreach(clusterroot) do n
            push!(mults, m[n])
        end
    end

    mults
end
