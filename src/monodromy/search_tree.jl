export SearchTree, iscontained, add!

const DEFAULT_CAPACITY = Ref(7)

struct SearchBlock{T}
    elements::Vector{Int32}
    children::Vector{Union{Nothing, SearchBlock{T}}}
    capacity::Int

    distances_cache::Vector{Tuple{T, Int}}
end

function SearchBlock(::Type{T}; capacity = DEFAULT_CAPACITY[]) where T
    children = Vector{Union{Nothing, SearchBlock{T}}}(nothing, capacity)
    distances_cache = Vector{Tuple{T, Int}}()
    SearchBlock(Int32[], children, capacity, distances_cache)
end

function SearchBlock(::Type{T}, index::Int; kwargs...) where T
    block = SearchBlock(T; kwargs...)
    push!(block.elements, index)
    block
end

function iscontained(block::SearchBlock, x::V, tol::Real, points::Vector{V}) where V
    if isempty(block.elements)
        return false
    end

    n = length(block.elements)
    d = distance(points[block.elements[1]], x)
    minᵢ = 1

    distances = block.distances_cache
    empty!(distances)
    push!(distances, (d, 1))

    if d < tol
        return true
    end

    for i ∈ 2:n
        dᵢ = distance(points[block.elements[i]], x)
        if dᵢ < d
            if dᵢ < tol
                return true
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
            if iscontained(block.children[minᵢ], x, tol, points)
                return true
            end
        else
            return false
        end
    end

    return false
end
iscontained(::Nothing, x::V, tol::Real, points::Vector{V}) where V = false

# This assumes that distances_cache is filled
function _insert!(block::SearchBlock{T}, index::Integer) where {T, V}
    if isempty(block.elements)
        push!(block.elements, index)
    end
    N = length(block.distances_cache)
    @assert N != 0

    dᵢ, minᵢ = block.distances_cache[end]
    # if not filled so far, just add it to the current block
    if length(block.elements) < block.capacity
        push!(block.elements, index)
    # we have no children so far, so create a new one
    elseif block.children[minᵢ] === nothing
        block.children[minᵢ] = SearchBlock(T, index; capacity=block.capacity)
    else # a block already exists, so recurse
        _insert!(block.children[minᵢ], index)
    end
end

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

function iscontained(tree::SearchTree{V}, x::V; tol::Real=1e-5) where V
    iscontained(tree.root, x, tol, tree.points)
end

function add!(tree::SearchTree{V}, x::V; tol::Real=1e-5) where V
    if iscontained(tree.root, x, tol, tree.points)
        return false
    end

    push!(tree.points, x)
    _insert!(tree.root, length(tree.points))
    true
end

function distance(x::AbstractVector{T}, y::AbstractVector{T}) where T
    n = length(x)
    @inbounds d = abs2(x[1] - y[1])
    @inbounds for i=2:n
        @fastmath d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
