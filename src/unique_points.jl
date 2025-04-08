export GroupActions,
    SymmetricGroup, UniquePoints, search_in_radius, add!, multiplicities, unique_points

################
# Group actions
################
"""
    GroupActions(actions::Function...)

Store a bunch of group actions `(f1, f2, f3, ...)`.
Each action has to return a tuple.
The actions are applied in the following sense
1) f1 is applied on the original solution `s`
2) f2 is applied on `s` and the results of 1
3) f3 is applied on `s` and the results of 1) and 2)
and so on

## Example
```julia-repl
julia> f1(s) = (s * s,);

julia> f2(s) = (2s, -s, 5s);

julia> f3(s) = (s + 1,);

julia> GroupActions(f1)(3)
(3, 9)

julia> GroupActions(f1, f2)(3)
(3, 9, 6, -3, 15, 18, -9, 45)

julia> GroupActions(f1,f2, f3)(3)
(3, 9, 6, -3, 15, 18, -9, 45, 4, 10, 7, -2, 16, 19, -8, 46)
```
"""
struct GroupActions{T<:Tuple}
    actions::T
end
GroupActions(::Nothing) = GroupActions(())
GroupActions(actions::GroupActions) = actions
GroupActions(actions::Function...) = GroupActions(actions)
GroupActions(actions) = GroupActions(actions...)

function (actions::GroupActions)(s)
    S = [s]
    T = typeof(s)
    apply_actions(actions, s) do sᵢ
        sⱼ = convert(T, sᵢ)
        if sⱼ != s
            push!(S, sⱼ)
        end
        false
    end
    S
end

apply_actions(cb, action::GroupActions, s) = _apply_actions(action.actions, s, cb)
@inline function _apply_actions(actions::Tuple, x, cb::F) where {F}
    f, rest = first(actions), Base.tail(actions)
    y = f(x)
    if isa(x, AbstractVector{<:Number}) && isa(y, AbstractVector{<:Number})
        cb(y) && return true
        if _apply_actions(rest, y, cb)
            return true
        end
    else
        for yᵢ in f(x)
            cb(yᵢ) && return true
            if _apply_actions(rest, yᵢ, cb)
                return true
            end
        end
    end
    _apply_actions(rest, x, cb)
end
@inline _apply_actions(::Tuple{}, s, cb) = false

apply_actions(cb::G, actions::F, s) where {G,F} = actions(cb, s)

# Implemented group actions

"""
    SymmetricGroup(n)

Group action of the symmetric group S(n).
"""
struct SymmetricGroup
    permutations::Vector{Vector{Int}}
end
SymmetricGroup(N::Int) = SymmetricGroup(permutations(N))
function permutations(N::Int)
    s = Vector(1:N)
    perms = [copy(s)]
    while true
        i = N - 1
        while i >= 1 && s[i] >= s[i+1]
            i -= 1
        end
        if i > 0
            j = N
            while j > i && s[i] >= s[j]
                j -= 1
            end
            s[i], s[j] = s[j], s[i]
            reverse!(s, i + 1)
        else
            s[1] = N + 1
        end

        s[1] > N && break
        push!(perms, copy(s))
    end
    perms
end

Base.eltype(::Type{SymmetricGroup}) = Vector{Int}
Base.length(p::SymmetricGroup) = length(p.permutations)
Base.iterate(p::SymmetricGroup) = iterate(p.permutations)
Base.iterate(p::SymmetricGroup, s) = iterate(p.permutations, s)


#############
# UniquePoints
#############

"""
    UniquePoints{T, Id, M}

A data structure for assessing quickly whether a point is close to an indexed point as
determined by the given distance function `M`. 
The indexed points are only stored by their identifiers `Id`. 
`triangle_inequality` should be set to `true`, if the distance function satisfies the triangle inequality. 
Otherwise, it should be set to `false`. If `triangle_inequality` is nothing the algorithm will try to detect whether the triangle is satisfied.

    UniquePoints(v::AbstractVector{T}, id::Id;
                    distance = InfNorm(),
                    triangle_inequality = nothing,
                    group_actions = nothing)

    

Initialize the data structure. This *does not* initialize the data structure with the point.


## Example

```julia
x = randn(ComplexF64, 4)
permutation(x) = ([x[2]; x[1]; x[3]; x[4]],)
group_actions = GroupActions(permutation)
X = group_actions(x)

# without group actions
unique_points = UniquePoints(x, 1)
HC.add!.(unique_points, X, 1:length(X), 1e-5)
length(unique_points) # 2

unique_points = UniquePoints(x, 1, group_actions = group_actions)
HC.add!.(unique_points, X, 1:length(X), 1e-5)
length(unique_points) # 1
```
"""
struct UniquePoints{T,Id,M,MaybeGA}
    tree::VoronoiTree{T,Id,M}
    group_actions::MaybeGA
    zero_vec::Vector{T}
end

function UniquePoints(
    v::AbstractVector,
    id;
    distance = InfNorm(),
    triangle_inequality = nothing,
    group_action = nothing,
    group_actions = isnothing(group_action) ? nothing : GroupActions(group_action),
)
    if (group_actions isa Tuple) || (group_actions isa AbstractVector)
        group_actions = GroupActions(group_actions)
    end

    d = distance

    if isnothing(triangle_inequality)
        if typeof(d) <: AbstractNorm
            triangle_inequality = true
        else
            n = length(v)
            v₁ = randn(ComplexF64, n)
            v₂ = randn(ComplexF64, n)
            v₃ = randn(ComplexF64, n)
            if d(v₁, v₂) ≤ d(v₁, v₃) + d(v₃, v₂) &&
               d(v₁, 4 .* v₁) ≤ d(v₁, 2 .* v₁) + d(2 .* v₁, 4 .* v₁)
                triangle_inequality = true
            else
                triangle_inequality = false
            end
        end
    end

    tree = VoronoiTree(v, id; distance = d, triangle_inequality = triangle_inequality)
    UniquePoints(tree, group_actions, zeros(eltype(v), length(v)))
end

function Base.show(io::IO, UP::UniquePoints)
    print(io, typeof(UP), " with ", length(UP.tree), " points")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::UniquePoints) = x
Base.length(UP::UniquePoints) = length(UP.tree)
Base.collect(UP::UniquePoints) = collect(UP.tree)
Base.broadcastable(UP::UniquePoints) = Ref(UP)
function Base.empty!(UP::UniquePoints{T,Id,M}) where {T,Id,M}
    empty!(UP.tree)
    UP
end

"""
    search_in_radius(unique_points, v, tol)

Search whether `unique_points` contains a point `p` with distances at most `tol` from `v`.
Returns `nothing` if no point exists, otherwise the identifier of `p` is returned.
"""
function search_in_radius(UP::UniquePoints{T,Id,M,GA}, v, tol::Real) where {T,Id,M,GA}
    id = search_in_radius(UP.tree, v, tol)
    if isnothing(id) && !isnothing(UP.group_actions)
        let actions = UP.group_actions::GA
            apply_actions(actions, v) do w
                id′ = search_in_radius(UP.tree, w, tol)
                if !isnothing(id′)
                    id = id′
                    return true
                end
                false
            end
        end
    end
end

"""
    add!(unique_points, v, id; atol = 1e-14, rtol = sqrt(eps()))
    add!(unique_points, v, id, atol)

Search whether `unique_points` contains a point `p` with distances at most
`max(atol, norm(v)rtol)` from `v`. If this is the case the identifier of `p` and `false` is
returned. Otherwise `(id, true)` is returned.
"""
function add!(UP::UniquePoints{T,Id,M,GA}, v, id::Id, tol::Real) where {T,Id,M,GA}
    found_id = search_in_radius(UP.tree, v, tol)
    if isnothing(found_id)
        if isnothing(UP.group_actions)
            insert!(UP.tree, v, id; use_distances = true)
            return (id, true)
        else
            let actions = UP.group_actions::GA
                apply_actions(actions, v) do w
                    found_id′ = search_in_radius(UP.tree, w, tol)
                    if !isnothing(found_id′)
                        found_id = found_id′
                        return true
                    end
                    false
                end
            end
            if isnothing(found_id)
                insert!(UP.tree, v, id)
                return (id, true)
            else
                return (found_id::Id, false)
            end
        end
    else
        return (found_id::Id, false)
    end
end
function add!(
    UP::UniquePoints{T,Id,M,GA},
    v,
    id::Id;
    atol::Float64 = 1e-14,
    rtol::Float64 = sqrt(eps()),
) where {T,Id,M,GA}
    n = UP.tree.distance(v, UP.zero_vec)
    rad = max(atol, rtol * n)
    add!(UP, v, id, rad)
end

####################
## Multiplicities ##
####################
"""
    multiplicities(vectors; distance = InfNorm(), atol = 1e-14, rtol = 1e-8, kwargs...)

Returns a `Vector{Vector{Int}}` `v`. Each vector `w` in 'v' contains all indices `i`,`j`
such that `w[i]` and `w[j]` have `distance` at most `max(atol, rtol * distance(0,w[i]))`.
The remaining `kwargs` are things that can be passed to [`UniquePoints`](@ref).

```julia-repl
julia> multiplicities([[1,0.5], [1,0.5], [1,1]])
[[1,2]]
```
This is the same as
```julia
multiplicities([[1,0.5], [1,0.5], [1,1]]; distance=(x,y) -> LinearAlgebra.norm(x-y))
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
function multiplicities(f::F, v; distance = InfNorm(), kwargs...) where {F<:Function}
    isempty(v) && return Vector{Vector{Int}}()
    _multiplicities(f, v, distance; kwargs...)
end
function _multiplicities(
    f::F,
    V,
    distance;
    atol::Float64 = 1e-14,
    rtol::Float64 = 1e-8,
    kwargs...,
) where {F<:Function}
    unique_points = UniquePoints(f(first(V)), 1; distance = distance, kwargs...)
    mults = Dict{Int,Vector{Int}}()
    for (i, vᵢ) in enumerate(V)
        wᵢ = f(vᵢ)
        k, new_point = add!(unique_points, wᵢ, i; atol = atol, rtol = rtol)
        if !new_point
            if haskey(mults, k)
                push!(mults[k], i)
            else
                mults[k] = [k, i]
            end
        end
    end
    collect(values(mults))
end
"""
    unique_points(vectors; distance = InfNorm(), atol = 1e-14, rtol = 1e-8, kwargs...)

Returns all elements in `vector` for which two elements have `distance` at most `max(atol, rtol * distance(0,w[i]))`.
Note that the output can depend on the order of elements in vectors.
The remaining `kwargs` are things that can be passed to [`UniquePoints`](@ref).
"""
function unique_points(V; distance = InfNorm(), atol = 1e-14, rtol = 1e-8, kwargs...)
    unique_points = UniquePoints(first(V), 1; distance = distance, kwargs...)
    out = Vector{eltype(V)}()
    for (i, vᵢ) in enumerate(V)
        _, new_point = add!(unique_points, vᵢ, i; atol = atol, rtol = rtol)
        if new_point
            push!(out, vᵢ)
        end
    end
    out
end
