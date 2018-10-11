"""
    Node(p::SVector, x::AbstractVector; store_points=true, is_main_node=false)

Create a node with a parameter from the same type as `p` and expecting
points with tthe same type as `x`. If `stores_points` is `true` a `UniquePoints`
data structure is allocated to keep track of all known solutions of this node.
"""
struct Node{N, T, UP<:UniquePoints}
    p::SVector{N, T} # maybe just a Vector?
    points::Union{Nothing, UP}
    # Metadata / configuration
    main_node::Bool
end

function Node(p::SVector{N, T}, x::AbstractVector; store_points=true, is_main_node=false) where {N, T}
    uniquepoints = UniquePoints(typeof(x))
    if store_points == false
        Node{N, T, typeof(uniquepoints)}(p, nothing, is_main_node)
    else
        Node(p, uniquepoints, is_main_node)
    end
end
function Node(p::SVector{N, T}, node::Node{N,T,UP}; store_points=true, is_main_node=false) where {N, T, UP}
    if store_points == false
        Node{N, T, UP}(p, nothing, is_main_node)
    else
        Node{N, T, UP}(p, UniquePoints(UP), is_main_node)
    end
end

apply_group_action(node::Node) = node.points !== nothing
add!(node::Node, x; kwargs...) = Utilities.add!(node.points, x; kwargs...)
function iscontained(node::Node, x; kwargs...)
    if node.points === nothing
        false
    else
        Utilities.iscontained(node.points, x; kwargs...)
    end
end
function unsafe_add!(node::Node, x; kwargs...)
    if node.points !== nothing
        Utilities.unsafe_add!(node.points, x; kwargs...)
    end
    nothing
end

struct Edge
    start::Int
    target::Int
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end
function Edge(start::Int, target::Int; usegamma=false)
    if usegamma
        γ = (randn(ComplexF64), randn(ComplexF64))
    else
        γ = nothing
    end
    Edge(start, target, γ)
end

function regenerate(edge::Edge)
    if edge.γ === nothing
        edge
    else
        γ = (randn(ComplexF64), randn(ComplexF64))
        Edge(edge.start, edge.target, γ)
    end
end


struct Graph{N<:Node}
    nodes::Vector{N}
    loop::Vector{Edge}
end

function Graph(p₁::SVector, x₁::AbstractVector, nnodes::Int, options::Options; usegamma=true)
    n₁ = Node(p₁, x₁, is_main_node = true)
    nodes = [n₁]
    store_points = options.group_action_on_all_nodes && has_group_actions(options)
    for i = 2:nnodes
        p = options.parameter_sampler(p₁)
        push!(nodes, Node(p, x₁, store_points=store_points))
    end

    loop = map(i -> Edge(i - 1, i; usegamma=usegamma), 2:nnodes)
    push!(loop, Edge(nnodes, 1; usegamma=usegamma))

    Graph(nodes, loop)
end


function add_initial_solutions!(G::Graph, solutions::Vector; kwargs...)
    for s ∈ solutions
        add!(G.nodes[1], s; kwargs...)
    end
    G
end

solutions(G::Graph) = G.nodes[1].points

mainnode(G::Graph) = G.nodes[1]
nodes(G::Graph) = G.nodes
loop(G::Graph) = G.loop

nextedge(G::Graph, edge::Edge) = G.loop[edge.target]

function regenerate!(G::Graph, parameter_sampler::Function, stats::Statistics)
    main = mainnode(G)

    # The first node is the main node and doesn't get touched
    for i ∈ 2:length(G.nodes)
        G.nodes[i] = Node(parameter_sampler(main.p), G.nodes[i])
    end
    G.loop .= regenerate.(G.loop)
    generated_parameters!(stats, length(main.points)) # bookkeeping
end

@inline function set_parameters!(tracker::PathTracking.PathTracker, e::Edge, G::Graph)
    H = Homotopies.basehomotopy(tracker.homotopy)
    if !(H isa Homotopies.ParameterHomotopy)
        error("Base homotopy is not a ParameterHomotopy")
    end
    p = (G.nodes[e.start].p, G.nodes[e.target].p)
    Homotopies.set_parameters!(H, p, e.γ)
end


function track(tracker, x::AbstractVector, edge::Edge, G::Graph, stats::Statistics)
    set_parameters!(tracker, edge, G)
    retcode = PathTracking.track!(tracker, x, 1.0, 0.0)
    pathtracked!(stats, retcode)
    iters = PathTracking.curriters(tracker)
    y = ProjectiveVectors.similar_affine(x, PathTracking.currx(tracker))
    return y, retcode
end
