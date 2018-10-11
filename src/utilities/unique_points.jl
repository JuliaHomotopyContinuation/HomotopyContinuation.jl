export UniquePoints, iscontained, add!, points

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

function iscontained(block::SearchBlock, x::AbstractVector, tol::Real, points::Vector, threadid=Threads.threadid())
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
iscontained(::Nothing, x::AbstractVector, tol::Real, points::Vector, threadid) = NOT_FOUND

# This assumes that distances_cache is filled
function _insert!(block::SearchBlock{T}, index::Integer, threadid=Threads.threadid()) where {T, V}
    if isempty(block.elements)
        push!(block.elements, index)
        return
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
    nothing
end

function Base.empty!(block::SearchBlock)
    empty!(block.elements)
    block.children .= nothing

    nothing
end

#############
# UniquePoints
#############

"""
    UniquePoints{V<:AbstractVector, T}

A data structure which holds points of type `V` where `T=real(eltype(V))`. This data structure
provides an efficient (poly)logarithmic check whether a point already exists where
two points are considered equal if their 2-norm is less than a provided tolerance `tol`.


    UniquePoints(v::AbstractVector{<:Number})

Initialize the data structure with just one data point `v`.


    UniquePoints(V::Vector{<:AbstractVector{<:Number}}; tol=1e-5)

Initialize the data structure with all points in `v`. These are added in order
by [`add!`](@ref) with the given tolerance `tol`.
"""
struct UniquePoints{V<:AbstractVector, T}
    root::SearchBlock{T}
    points::Vector{V}
end

UniquePoints(v::Type{<:UniquePoints{V}}) where V = UniquePoints(V)
function UniquePoints(::Type{V}) where {T<:Number, V<:AbstractVector{T}}
    root = SearchBlock(real(T))
    points = Vector{V}()
    UniquePoints(root, points)
end
function UniquePoints(v::AbstractVector{T}) where {T<:Number}
    root = SearchBlock(real(T), 1)
    points = [v]
    UniquePoints(root, points)
end

function UniquePoints(v::AbstractVector{<:AbstractVector}; kwargs...)
    data = UniquePoints(v[1])
    for i = 2:length(v)
        add!(data, v[i]; kwargs...)
    end
    data
end

function Base.similar(data::UniquePoints{V, T}) where {V, T}
    root = SearchBlock(T)
    points = Vector{V}()
    UniquePoints(root, points)
end

"""
    points(data::UniquePoints)

Return the points stored in `data`.
"""
points(data::UniquePoints) = data.points

Base.getindex(data::UniquePoints, i::Integer) = data.points[i]


"""
    iscontained(data::UniquePoints{V}, x::V; tol=1e-5)::Bool

Check whether `x` is contained in the `data` by using the tolerance `tol` to decide for duplicates.

    iscontained(data::UniquePoints{V}, x::V, Val(true); tol=1e-5)::Int

If `x` is contained in `data` by using the tolerance `tol` return the index
of the data point which already exists. If the data point is not existing `-1`
is returned.
"""
function iscontained(data::UniquePoints, x::AbstractVector, ::Val{Index}=Val{false}(); tol::Real=1e-5) where {Index}
    if Index
        iscontained(data.root, x, tol, data.points)
    else
        iscontained(data.root, x, tol, data.points) ≠ NOT_FOUND
    end
end

"""
    add!(data::UniquePoints{V}, x::V; tol=1e-5)::Bool

Add `x` to `data` if it doesn't already exists by using the tolerance `tol` to decide for duplicates.

    add!(data::UniquePoints{V}, x::V, Val(true); tol=1e-5)::Int

If `x` is contained in `data` by using the tolerance `tol` to decide for duplicates return the index
of the data point which already exists. If the data point is not existing add it to `x` and
return `-1`. The element will be the last element of `points(data)`.
"""
function add!(data::UniquePoints, x::AbstractVector, ::Val{Index}=Val{false}(); tol::Real=1e-5) where {Index}
    if Index
        idx = iscontained(data.root, x, tol, data.points)
        if idx ≠ NOT_FOUND
            return idx
        end
        push!(data.points, x)
        _insert!(data.root, length(data.points))
        NOT_FOUND
    else
        if iscontained(data.root, x, tol, data.points) ≠ NOT_FOUND
            return false
        end

        push!(data.points, x)
        _insert!(data.root, length(data.points))
        true
    end
end

"""
    unsafe_add!(data::UniquePoints{V}, x::V)::Bool

Similarly to [`add!`](@ref) but assumes that it was already checked that there is no
duplicate with [`iscontained`](@ref). *This has to be called directly after `iscontained`
with the same value of `x`*.
"""
function unsafe_add!(data::UniquePoints, x::AbstractVector) where {Index}
    push!(data.points, x)
    _insert!(data.root, length(data.points))

end

Base.length(data::UniquePoints) = length(data.points)

Base.iterate(data::UniquePoints) = iterate(data.points)
Base.iterate(data::UniquePoints, state) = iterate(data.points, state)
Base.eltype(::Type{<:UniquePoints{V}}) where {V} = V
Base.size(data::UniquePoints) = size(data.points)

"""
    empty!(data::UniquePoints)

Remove all points from `data`.
"""
function Base.empty!(data::UniquePoints)
    empty!(data.root)
    empty!(data.points)
end

function distance(x::AbstractVector{T}, y::AbstractVector{T}) where T
    n = length(x)
    @inbounds d = abs2(x[1] - y[1])
    @inbounds for i=2:n
        @fastmath d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
