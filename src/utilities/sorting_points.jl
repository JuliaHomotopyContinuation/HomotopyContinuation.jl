export UniquePoints, multiplicities, iscontained, add!, simple_add!, empty!, points, unique_points

const DEFAULT_CAPACITY = Ref(7) # Determined by testing a couple of different values
const NOT_FOUND = -1
const NOT_FOUND_AND_REAL = -2

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

function iscontained(block::SearchBlock{T}, x::NTuple{N,S}, tol::Real, points::Vector, distance::F, threadid=Threads.threadid()) where {T, N, S, F<:Function}
    iscontained(block, SVector{N,S}(x), tol, points, distance)
end

function iscontained(block::SearchBlock{T}, x::AbstractVector{S}, tol::Real, points::Vector, distance::F, threadid=Threads.threadid()) where {T, S, F<:Function}
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
    retidx = iscontained(block.children[m₁[2]], x, tol, points, distance, threadid)
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
    retidx = iscontained(block.children[m₂[2]], x, tol, points, distance, threadid)
    if retidx ≠ NOT_FOUND
        distances[1] = m₁ # we rely on the distances for look up, so place at the first place the smallest element
        return retidx
    end

    if m₃[1] - m₁[1] > 2tol
        distances[1] = m₁
        return NOT_FOUND # Checked first and second case
    end

    # Since we know als the third element, let's check it
    retidx = iscontained(block.children[m₃[2]], x, tol, points, distance, threadid)
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
            retidx = iscontained(block.children[i], x, tol, points, distance, threadid)
            if retidx ≠ NOT_FOUND
                return retidx
            end
        else
            break
        end
    end

    return NOT_FOUND
end
iscontained(::Nothing, x::AbstractVector, tol::Real, points::Vector, distance, threadid) = NOT_FOUND

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
    UniquePoints{V<:AbstractVector, T, F<:Function}

A data structure which holds points of type `V` where `T=real(eltype(V))`. This data structure
provides an efficient (poly)logarithmic check whether a point already exists where
two points `u,v` are considered equal if `F(u,v)<tol`, where `tol` is a tolerance provided through the [`add!`](@ref) function.


    UniquePoints(v::AbstractVector{<:Number}, distance::F)

Initialize the data structure with just one data point `v`.


    UniquePoints(V::Vector{<:AbstractVector{<:Number}}, distance::F; tol=1e-5)

Initialize the data structure with all points in `v`. These are added in order
by [`add!`](@ref) with the given tolerance `tol`. In particular, 'UniquePoints' structure will contain only points for which the pairwise distance given by `F` is less than `tol`.

    UniquePoints(v; kwargs...) = UniquePoints(v, euclidean_distance; kwargs...)

If `F` is not specialized, [`euclidean_distance`](@ref) is used.

Optional keywords:

* `check_real=true` adds real from points from group orbits (if they exist). The default is `check_real=true`.
* The user can use `group_action=foo` or, if there is more than one group acting, `group_actions=[foo, bar]`. Then, points that are in the same group orbit are considered equal. See [`GroupActions`](@ref) for details regarding the application rules.

## Examples
```julia-repl
julia> points(UniquePoints([[1.0,0.5], [1.0,0.5], [0.5,1.0]]))
2-element Array{Array{Float64,1},1}:
 [1.0, 0.5]
 [0.5, 1.0]

julia> points(UniquePoints([[1.0,0.5], [1.0,0.5], [0.5,1.0]], group_action = x -> [x[2],x[1]]))
1-element Array{Array{Float64,1},1}:
 [1.0, 0.5]
```
"""
struct UniquePoints{V<:AbstractVector, T, F<:Function, MaybeGA<:Union{Nothing, GroupActions}}
    root::SearchBlock{T}
    points::Vector{V}
    distance_function::F
    group_actions::MaybeGA
    check_real::Bool
end

function UniquePoints(::Type{V}, distance::F;
    group_action=nothing,
    group_actions=group_action === nothing ? nothing : GroupActions(group_action),
    check_real::Bool=true) where {T<:Number, V<:AbstractVector{T}, F<:Function}

    if group_actions isa GroupActions
       actions = group_actions
   elseif group_actions === nothing
       actions = nothing
    else
       actions = GroupActions(group_actions)
    end

    root = SearchBlock(real(T))
    points = Vector{V}()
    UniquePoints(root, points, distance, actions, check_real)
end

function UniquePoints(v::AbstractVector{<:Number}, distance::Function; kwargs...)
    data = UniquePoints(typeof(v), distance; kwargs...)
    add!(data, v)
    data
end

function UniquePoints(v::AbstractVector{<:AbstractVector}, distance::F; tol::Real=1e-5, kwargs...) where {F<:Function}
    data = UniquePoints(eltype(v), distance; kwargs...)
    add!(data, v; tol=tol)
end
UniquePoints(v::Type{<:UniquePoints{V}}, distance::F; kwargs...) where {V, F<:Function} = UniquePoints(V, distance; kwargs...)
UniquePoints(v; distance=euclidean_distance, kwargs...) = UniquePoints(v, distance; kwargs...)

function Base.show(io::IO, data::UniquePoints)
    print(io, typeof(data), " with ", length(points(data)), " points")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::UniquePoints) = x

"""
    points(data::UniquePoints)

Return the points stored in `data`.
"""
points(data::UniquePoints) = data.points


Base.getindex(data::UniquePoints, i::Integer) = data.points[i]


"""
    iscontained(data::UniquePoints{V}, x::V; tol=1e-5)::Bool

Check whether `x` is contained in the `data` by using the tolerance `tol` to decide for duplicates.

    iscontained(data::UniquePoints{V}, x::V, Val{true}(); tol=1e-5)::Int

If `x` is contained in `data` by using the tolerance `tol` return the index
of the data point which already exists. If the data point is not existing `-1`
is returned. If `data` has the option `check_real` enabled, a `-2` will be returned once a real vector was added.
"""
function iscontained(data::UniquePoints, x, val=Val{false}(); tol::Real=1e-5)
    iscontained(data, x, val, tol)
end
function iscontained(data::UniquePoints, x::NTuple{N,T}, val, tol::Real) where {N, T}
    iscontained(data, SVector{N,T}(x), val, tol)
end
function iscontained(data::UniquePoints{T}, x::AbstractVector{S}, ::Val{Index}, tol::Real) where {T,S,Index}
    index = iscontained(data.root, x, tol, data.points, data.distance_function)
    if index == NOT_FOUND
        if data.group_actions !== nothing # extra if statement since inference cannot look through &&
            apply_actions(data.group_actions, x) do y
                k = iscontained(data.root, y, tol, data.points, data.distance_function)
                if k ≠ NOT_FOUND
                    index = k
                    return true
                end
                false
            end
        end
    end

    if Index
        return index
    else
        return index ≠ NOT_FOUND
    end
end

"""
    add!(data::UniquePoints{V}, x::V; tol=1e-5)::Bool

Add `x` to `data` if it doesn't already exists by using the tolerance `tol` to decide for duplicates.

    add!(data::UniquePoints{V}, x::V, Val(true); tol=1e-5)::Int

If `x` is contained in `data` by using the tolerance `tol` to decide for duplicates return the index
of the data point which already exists. If the data point is not existing add it to `data` and
return `-1`. If `data` has the option `check_real` enabled, a `-2` will be returned once a real vector was added. The element will be the last element of `points(data)`.
"""
function add!(data::UniquePoints, x::AbstractVector{<:Number}, ::Val{true}; tol::Real=1e-5)
    idx = iscontained(data, x, Val(true), tol)
    idx == NOT_FOUND || return idx

    if !data.check_real
        simple_add!(data, x, tol)
        return NOT_FOUND
    else
        if isrealvector(x)
            simple_add!(data, x, tol)
            return NOT_FOUND_AND_REAL
        elseif data.group_actions !== nothing
            not_found_and_real = false
            apply_actions(data.group_actions, x) do y
                if isrealvector(y)
                    simple_add!(data, y, tol)
                    not_found_and_real = true
                    return true
                end
                false
            end
            not_found_and_real && return NOT_FOUND_AND_REAL
        end
        simple_add!(data, x, tol)
        return NOT_FOUND
    end
end
function add!(data::UniquePoints, x::AbstractVector{<:Number}, ::Val{false}=Val(false); tol::Real=1e-5)
    idx = add!(data, x, Val(true); tol=tol)
    idx == NOT_FOUND || idx == NOT_FOUND_AND_REAL
end
function add!(data::UniquePoints, v::AbstractVector{<:AbstractVector}, val::Val=Val(false); tol::Real=1e-5)
    for vᵢ in v
        add!(data, vᵢ; tol=tol)
    end
    data
end

"""
    simple_add!(data::UniquePoints{V}, x::V, tol::Real)::Bool

Similarly to [`add!`](@ref) but does not apply any group actions.
If the data point is not existing add it to `data` and return `-1`. Otherwise
the index of `x` in `data.points` is returned.
"""
function simple_add!(data::UniquePoints, x::AbstractVector, tol::Real)
    idx = iscontained(data.root, x, tol, data.points, data.distance_function)
    if idx ≠ NOT_FOUND
        return idx
    end
    push!(data.points, x)
    _insert!(data.root, length(data.points))
    NOT_FOUND
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



"""
    multiplicities(vectors; distance=euclidean_distance, tol::Real = 1e-5, kwargs...)

Returns an array of arrays of integers. Each vector `w` in 'v' contains all indices `i,j` such that `w[i]` and `w[j]` have `distance` at most tol.

Optional keywords:

* `check_real=true` adds real from points from group orbits (if they exist) to the [`UniquePoints`](@ref) data structure used internally. The default is `check_real=false`.
* The user can use `group_action=foo` or, if there is more than one group acting, `group_actions=[foo, bar]`. Then, points that are in the same group orbit are considered equal. See [`GroupActions`](@ref) for details regarding the application rules.


```julia-repl
julia> multiplicities([[1,0.5]; [1,0.5]; [1,1]])
[[1,2]]
```
This is the same as
```julia
multiplicities([[1,0.5]; [1,0.5]; [1,1]]; distance=(x,y) -> LinearAlgebra.norm(x-y))
```
Here is an example for using group actions.
```julia-repl
julia> X = [[1, 2, 3, 4], [2,1,3,4], [1,2,4,3], [2,1,4,3]]
julia> permutation(x) = [x[2], x[1], x[3], x[4]]
julia> m = multiplicities(X, group_action = permutation)
[[1,2], [3,4]]
```
"""
multiplicities(v; kwargs...) = multiplicities(identity, v; kwargs...)
function multiplicities(f, v; distance=default_distance(f,v), kwargs...) where {F<:Function}
    _multiplicities(f, v, distance; kwargs...)
end
default_distance(f, v) = isa(f(first(v)), PVector) ? fubini_study : euclidean_distance
function _multiplicities(f, v, distance::F; tol::Float64=1e-5, check_real=false, kwargs...) where {F<:Function}
    mults = Dict{Int, Vector{Int}}()
    positions = Vector{Int32}()
    x₀ = f(first(v))
    T = typeof(similar(x₀, promote_type(eltype(x₀), Float64)))
    data = UniquePoints(T, distance; check_real=check_real, kwargs...)
    for (i, vᵢ) in enumerate(v)
        k = add!(data, f(vᵢ), Val(true); tol = tol)
        if k > 0
            if haskey(mults, positions[k])
                push!(mults[positions[k]], i)
            else
                mults[positions[k]] = [positions[k], i]
            end
        else
            push!(positions, i)
        end
    end
    collect(values(mults))
end

"""
    unique_points(points::AbstractVector{<:AbstractVector}; options...)

Compute all unique points with respect to the given options. See [`UniquePoints`](@ref)
for possible options. In particular, it is possible to pass group actions.

## Example
```julia-repl
julia> unique_points([[1.0,0.5], [1.0,0.5], [0.5,1.0]])
2-element Array{Array{Float64,1},1}:
 [1.0, 0.5]
 [0.5, 1.0]

julia> unique_points([[1.0,0.5], [1.0,0.5], [0.5,1.0]]; group_action = x -> [x[2],x[1]])
1-element Array{Array{Float64,1},1}:
 [1.0, 0.5]
"""
unique_points(v::AbstractVector{<:AbstractVector}; kwargs...) = points(UniquePoints(v; kwargs...))
