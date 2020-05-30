mutable struct VTNode{T,Id}
    nentries::Int
    values::Matrix{T}
    ids::Vector{Id}
    children::Vector{VTNode{T,Id}}
    distances::Vector{Tuple{Float64,Int}}
end

function VTNode(::Type{T}, ::Type{Id}, d::Int; capacity::Int) where {T,Id}
    values = Matrix{T}(undef, d, capacity)
    ids = Vector{Id}(undef, capacity)
    children = Vector{VTNode{T,Id}}(undef, capacity)
    distances = Vector{Tuple{Float64,Int}}(undef, capacity)
    VTNode(0, values, ids, children, distances)
end

function VTNode(v::AbstractVector{T}, id::Id; kwargs...) where {T,Id}
    node = VTNode(T, Id, length(v); kwargs...)
    node.nentries = 1
    node.values[:, 1] .= v
    node.ids[1] = id
    node
end

Base.length(node::VTNode) = node.nentries
Base.isempty(node::VTNode) = length(node) == 0
capacity(node::VTNode) = length(node.children)

function compute_distances!(node, x, metric::M) where {M}
    for j = 1:length(node)
        node.distances[j] = (metric(x, view(node.values, :, j)), j)
    end
    node.distances
end

function search_in_radius(
    node::VTNode{T,Id},
    x::AbstractVector,
    tol::Real,
    metric::M,
) where {T,Id,M}
    !isempty(node) || return nothing

    n = length(node)

    # We compute now the distance to every other element in the node
    # If the distance to one element is smaller than tol, we are done.
    # Otherwise, we need to look into those children whose distance
    # is closest (let's say with distance `d`) and all those
    # children for whose distance `dᵢ` we have |d-dᵢ| < 2*tol
    # This fact can be shown by using the triangle inequality.
    # We therefore simply compute the distances to every point and sort
    # afterwards

    # We go through the elements and compute the distances,
    # keeping track of the three smallest elements
    m₁ = m₂ = m₃ = (Inf, 1)
    # we have a distances cache per thread
    distances = compute_distances!(node, x, metric)
    for i = 1:n
        dᵢ = first(distances[i])
        # early exit
        if dᵢ < tol
            # we rely on the distances for insertion, so place at the first place the smallest element
            distances[1] = (dᵢ, i)
            return node.ids[i]
        end

        # check three smallest elements and update if necessary
        if dᵢ < first(m₁)
            m₃ = m₂
            m₂ = m₁
            m₁ = (dᵢ, i)
        elseif dᵢ < first(m₂)
            m₃ = m₂
            m₂ = (dᵢ, i)
        elseif dᵢ < first(m₃)
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
    if isassigned(node.children, last(m₁))
        retid = search_in_radius(node.children[last(m₁)], x, tol, metric)
        if !isnothing(retid)
            # we rely on the distances for insertion, so place the smallest element first
            distances[1] = m₁
            return retid::Id
        end
    end
    # Case 1)
    if m₂[1] - m₁[1] > 2tol
        distances[1] = m₁
        return nothing # already checked first tree
    end
    # Case 2) We know m₂[1] - m₁[1] ≤ 2tol
    if isassigned(node.children, last(m₂))
        retid = search_in_radius(node.children[last(m₂)], x, tol, metric)
        if !isnothing(retid)
            # we rely on the distances for insertion, so place the smallest element first
            distances[1] = m₁
            return retid::Id
        end
    end

    if m₃[1] - m₁[1] > 2tol
        distances[1] = m₁
        return nothing # Checked first and second case
    end

    # Since we know also the third element, let's check it
    if isassigned(node.children, last(m₃))
        retid = search_in_radius(node.children[last(m₃)], x, tol, metric)
        if !isnothing(retid)
            # we rely on the distances for insertion, so place at the first place the smallest element
            distances[1] = m₁
            return retid::Id
        end
    end

    # Case 3)
    # We need to sort distances
    sort!(view(distances, 1:n), Base.Sort.InsertionSort, Base.Sort.By(first))

    # We can start at 4 since we already checked the smallest 3
    d = m₃[1]
    for k ∈ 4:n
        dᵢ, i = distances[k]
        if dᵢ - m₁[1] < 2tol
            if isassigned(node.children, i)
                retid = search_in_radius(node.children[i], x, tol, metric)
                if !isnothing(retid)
                    return retid::Id
                end
            end
        else
            break
        end
    end

    return nothing
end


function _insert!(
    node::VTNode{T,Id},
    v::AbstractVector,
    id::Id,
    metric::M;
    use_distances::Bool = false,
) where {T,Id,M}
    # if not filled so far, just add it to the current node
    if length(node) < capacity(node)
        k = (node.nentries += 1)
        node.values[:, k] .= v
        node.ids[k] = id
        return nothing
    end

    if use_distances
        dᵢ, minᵢ = first(node.distances)
    else
        compute_distances!(node, v, metric)
        dᵢ, minᵢ = findmin(node.distances)
    end

    if !isassigned(node.children, minᵢ)
        node.children[minᵢ] = VTNode(v, id; capacity = capacity(node))
    else # a node already exists, so recurse
        _insert!(node.children[minᵢ], v, id, metric)
    end

    nothing
end

function identifiers!(ids, node::VTNode)
    for i = 1:length(node)
        push!(ids, node.ids[i])
        if isassigned(node.children, i)
            identifiers!(ids, node.children[i])
        end
    end
    ids
end

"""
     VoronoiTree(
    v::AbstractVector{T},
    id::Id;
    metric = EuclideanNorm(),
    capacity = 8,
)

Construct a Voronoi tree data structure for vector `v` of element type `T` and with identifiers
`Id`. Each node has the given `capacity` and distances are measured by the given `metric`.
"""
mutable struct VoronoiTree{T,Id,M}
    root::VTNode{T,Id}
    nentries::Int
    metric::M
end

function VoronoiTree{T,Id}(d::Int; metric = EuclideanNorm(), capacity::Int = 8) where {T,Id}
    root = VTNode(T, Id, d; capacity = capacity)
    VoronoiTree(root, 0, metric)
end

function VoronoiTree(v::AbstractVector{T}, id::Id; kwargs...) where {T,Id}
    VoronoiTree{T,Id}(length(v); kwargs...)
end

Base.length(T::VoronoiTree) = T.nentries
Base.broadcastable(T::VoronoiTree) = Ref(T)

"""
    insert!(tree::VoronoiTree, v::AbstractVector, id)

Insert in the tree the point `v` with identifier `id`.
"""
function Base.insert!(
    tree::VoronoiTree{T,Id},
    v::AbstractVector,
    id::Id;
    use_distances::Bool = false,
) where {T,Id}
    _insert!(tree.root, v, id, tree.metric; use_distances = use_distances)
    tree.nentries += 1
    tree
end

"""
    search_in_radius(tree::VoronoiTree, v::AbstractVector, tol::Real)

Search whether the given `tree` contains a point `p` with distances at most `tol` from `v`.
Returns `nothing` if no point exists, otherwise the identifier of `p` is returned.
"""
function search_in_radius(tree::VoronoiTree, v::AbstractVector, tol::Real)
    search_in_radius(tree.root, v, tol, tree.metric)
end


"""
    add!(tree::VoronoiTree, v, id, tol)

Insert in the tree the point `v` with identifier `id` if `search_in_radius(tree v, tol)` is
`nothing`. If this is the case the tuple `(id, true)` is returned, otherwise the obtained
identifier and `false` is returned.
"""
function add!(tree::VoronoiTree{T,Id}, v, id::Id, tol::Real) where {T,Id}
    found_id = search_in_radius(tree, v, tol)
    if isnothing(found_id)
        insert!(tree, v, id; use_distances = true)
        return (id, true)
    else
        return (found_id::Id, false)
    end
end

function Base.collect(tree::VoronoiTree{T,Id}) where {T,Id}
    ids = Id[]
    identifiers!(ids, tree.root)
    ids
end
