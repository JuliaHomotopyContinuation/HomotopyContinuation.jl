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
    distances_cache = [distance_cache(T, capacity) for _=1:Threads.nthreads()]
    SearchBlock(Int32[], children, capacity, distances_cache)
end

function distance_cache(::Type{T}, capacity) where T
    [(typemax(T), i) for i=1:capacity]
end

function SearchBlock(::Type{T}, index::Int; kwargs...) where T
    block = SearchBlock(T; kwargs...)
    push!(block.elements, index)
    block
end

function iscontained(block::SearchBlock{T}, x::AbstractVector, tol::Real, points::Vector, threadid=Threads.threadid()) where T
    if isempty(block.elements)
        return NOT_FOUND
    end

    n = length(block.elements)

    # We compute now the distance to every other element in the block
    # If the distance to one element is smaller than tol, we are done.
    # Otherwise, we need to look into those child whose distance
    # is closest (let's say with distance `d`) and all those
    # children for whose distance `dᵢ` we have |d-dᵢ| < 2*tol
    # This fact can be shown by using the triangle inequality.
    # We therefore simply compute the distances to every point and sort
    # afterwards


    # We go through the elements and compute the distances,
    # keeping track of the three smallest elements
    m₁ = m₂ = m₃ = (typemax(T), 1)
    # we have a distances cache per thread
    distances = block.distances_cache[threadid]
    for i ∈ 1:n
        dᵢ = distance(points[block.elements[i]], x)
        # early exit
        if dᵢ < tol
             # we rely on the distances for look up, so place at the first place the smallest element
            distances[1] = (dᵢ, i)
            return block.elements[i]
        end

        distances[i] = (dᵢ, i)
        # check three smallest elements and update if necessary
        if dᵢ < m₁[1]
            m₃ = m₂
            m₂ = m₁
            m₁ = (dᵢ, i)
        elseif dᵢ < m₂[1]
            m₃ = m₂
            m₂ = (dᵢ, i)
        elseif dᵢ < m₃[1]
            m₃ = (dᵢ, i)
        end
    end

    # Now we computed all distances and also kept track of the two smallest elements
    # Now there are three cases
    # 1) m₁[1] + 2tol < m₂[1] -- > The element can only be in one subtree
    # 2) m₁[1] + 2tol < m₃[1] -- > The element can only be in the first or second subtree
    # 3) else -> The element can also be in more subtrees,
    #            we need to sort the distances vector and look through everything

    # Check smallest element first
    retidx = iscontained(block.children[m₁[2]], x, tol, points, threadid)
    if retidx ≠ NOT_FOUND
        distances[1] = m₁ # we rely on the distances for look up, so place at the first place the smallest element
        return retidx
    end
    # Case 1)
    if m₂[1] - m₁[1] > 2tol
        distances[1] = m₁
        return NOT_FOUND # already checked first tree
    end
    # Case 2) We know m₂[1] - m₁[1] ≤ 2tol
    retidx = iscontained(block.children[m₂[2]], x, tol, points, threadid)
    if retidx ≠ NOT_FOUND
        distances[1] = m₁ # we rely on the distances for look up, so place at the first place the smallest element
        return retidx
    end

    if m₃[1] - m₁[1] > 2tol
        distances[1] = m₁
        return NOT_FOUND # Checked first and second case
    end

    # Since we know als the third element, let's check it
    retidx = iscontained(block.children[m₃[2]], x, tol, points, threadid)
    if retidx ≠ NOT_FOUND
        distances[1] = m₁ # we rely on the distances for look up, so place at the first place the smallest element
        return retidx
    end

    # Case 3)
    # We need to sort distances
    sort!(distances, Base.Sort.InsertionSort, Base.Sort.By(first))

    # We can start at 4 since we already checked the smallest 3
    d = m₃[1]
    for k ∈ 4:n
        dᵢ, i = distances[k]
        if dᵢ - m₁[1] < 2tol
            retidx = iscontained(block.children[i], x, tol, points, threadid)
            if retidx ≠ NOT_FOUND
                return retidx
            end
        else
            break
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

    dᵢ, minᵢ = block.distances_cache[threadid][1]
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
