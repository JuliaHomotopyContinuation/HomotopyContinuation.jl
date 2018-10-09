export SearchTree, iscontained, add!

const DEFAULT_CAPACITY = Ref(7) # Determined by testing a couple of different values
const NOT_FOUND = -1

#############
# SearchBlock
#############
struct SearchBlock{T}
    elements::Vector{Int32}
    children::Vector{Union{Nothing, SearchBlock{T}}}
    capacity::Int
    distances_cache::Vector{Vector{Tuple{T, Int}}} # one per thread
end

function SearchBlock(::Type{T}; capacity = DEFAULT_CAPACITY[]) where T
    children = Vector{Union{Nothing, SearchBlock{T}}}(nothing, capacity)
    distances_cache = [Vector{Tuple{T, Int}}() for _=1:Threads.nthreads()]
    SearchBlock(Int32[], children, capacity, distances_cache)
end

function SearchBlock(::Type{T}, index::Int; kwargs...) where T
    block = SearchBlock(T; kwargs...)
    push!(block.elements, index)
    block
end

function iscontained(block::SearchBlock, x::V, tol::Real, points::Vector{V}, threadid=Threads.threadid()) where V
    if isempty(block.elements)
        return NOT_FOUND
    end

    n = length(block.elements)
    d = distance(points[block.elements[1]], x)
    minᵢ = 1

    distances = block.distances_cache[threadid]
    empty!(distances)
    push!(distances, (d, 1))

    if d < tol
        return block.elements[1]
    end

    for i ∈ 2:n
        dᵢ = distance(points[block.elements[i]], x)
        if dᵢ < d
            if dᵢ < tol
                return block.elements[i]
            end
            d = dᵢ
            minᵢ = i
            push!(distances, (dᵢ, minᵢ))
        end
    end

    N = length(distances)
    for i ∈ N:-1:1
        dᵢ, minᵢ = distances[i]
        if abs(d - dᵢ) < 2*tol # What is the correct constant?
            retidx = iscontained(block.children[minᵢ], x, tol, points, threadid)
            if retidx ≠ NOT_FOUND
                return retidx
            end
        else
            return NOT_FOUND
        end
    end

    return NOT_FOUND
end
iscontained(::Nothing, x::V, tol::Real, points::Vector{V}, threadid) where V = NOT_FOUND

# This assumes that distances_cache is filled
function _insert!(block::SearchBlock{T}, index::Integer, threadid=Threads.threadid()) where {T, V}
    if isempty(block.elements)
        push!(block.elements, index)
    end
    distances = block.distances_cache[threadid]
    N = length(distances)
    @assert N != 0

    dᵢ, minᵢ = distances[end]
    # if not filled so far, just add it to the current block
    if length(block.elements) < block.capacity
        push!(block.elements, index)
    # we have no children so far, so create a new one
    elseif block.children[minᵢ] === nothing
        block.children[minᵢ] = SearchBlock(T, index; capacity=block.capacity)
    else # a block already exists, so recurse
        _insert!(block.children[minᵢ], index, threadid)
    end
end

function Base.empty!(block::SearchBlock)
    empty!(block.elements)
    block.children .= nothing

    nothing
end

#############
# SearchTree
#############

"""
    SearchTree{V<:AbstractVector, T}

A data structure which holds all known solutions of type `V`. This data structure
provides an efficient (poly)logarithmic check whether a solution already exists where
two solutions are considered equal if their 2-norm is less than `tol`.
"""
struct SearchTree{V<:AbstractVector, T}
    root::SearchBlock{T}
    points::Vector{V}
end

function SearchTree(v::AbstractVector{T}) where {T<:Number}
    root = SearchBlock(real(T), 1)
    points = [v]
    SearchTree(root, points)
end

function SearchTree(v::AbstractVector{<:AbstractVector}; kwargs...)
    tree = SearchTree(v[1])
    for i = 2:length(v)
        add!(tree, v[i]; kwargs...)
    end
    tree
end

function iscontained(tree::SearchTree{V}, x::V, ::Val{Index}=Val{false}(); tol::Real=1e-5) where {V, Index}
    if Index
        iscontained(tree.root, x, tol, tree.points)
    else
        iscontained(tree.root, x, tol, tree.points) ≠ NOT_FOUND
    end
end

function add!(tree::SearchTree{V}, x::V, ::Val{Index}=Val{false}(); tol::Real=1e-5) where {V, Index}
    if Index
        idx = iscontained(tree.root, x, tol, tree.points)
        if idx ≠ NOT_FOUND
            return idx
        end
        push!(tree.points, x)
        _insert!(tree.root, length(tree.points))
        NOT_FOUND
    else
        if iscontained(tree.root, x, tol, tree.points) ≠ NOT_FOUND
            return false
        end

        push!(tree.points, x)
        _insert!(tree.root, length(tree.points))
        true
    end
end

Base.length(tree::SearchTree) = length(tree.points)

function Base.empty!(tree::SearchTree)
    empty!(tree.root)
    empty!(tree.points)
end

function distance(x::AbstractVector{T}, y::AbstractVector{T}) where T
    n = length(x)
    @inbounds d = abs2(x[1] - y[1])
    @inbounds for i=2:n
        @fastmath d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
