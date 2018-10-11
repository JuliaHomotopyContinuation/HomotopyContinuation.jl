"""
    Node(p₀::SVector, x₀::AbstractVector)

Create a node with a parameter from the same type as `p₀` and expecting
points with tthe same type as `x₀`.
"""
struct Node{N, T, UP<:UniquePoints}
    p::SVector{N, T} # maybe just a Vector?
    points::Union{Nothing, UP}
    # Metadata / configuration
    apply_group_action::Bool
    main_node::Bool
end

function Node(p::SVector{NParams, T}, x₀::AbstractVector; apply_group_action=true, main_node=false) where {NParams, T}
    Node(p, UniquePoints(typeof(x₀)), apply_group_action, main_node)
end

function regenerate(node::Node{N, T}, p; apply_group_action=node.apply_group_action) where {N, T}
    Node(p, similar(node.points), apply_group_action, false)
end

iscontained(node::Node, x; kwargs...) = iscontained(node.points, x; kwargs...)
add!(node::Node, x; kwargs...) = Utilities.add!(node.points, x; kwargs...)

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

function Graph(p₁::SVector, x₁::AbstractVector, nnodes::Int, options::Options; usegamma=true, apply_group_action_on=:main_node)
    main_group = tail_group = false
    if apply_group_action_on == :all_nodes
        main_group = tail_group = true
    elseif apply_group_action_on == :main_node
        main_group = true
    end

    n₁ = Node(p₁, x₁, apply_group_action = main_group, main_node = true)
    nodes = [n₁]
    for i = 2:nnodes
        p = options.parameter_sampler(p₁)
        push!(nodes, Node(p, x₁, apply_group_action=tail_group))
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


nodes(G::Graph) = G.nodes
loop(G::Graph) = G.loop

nextedge(G::Graph, edge::Edge) = G.loop[edge.target]

function regenerate!(G::Graph, parameter_sampler::Function)
    # The first node is fixed
    for i ∈ 2:length(G.nodes)
        G.nodes[i] = regenerate(G.nodes[i], parameter_sampler(G.nodes[1].p))
    end
    G.loop .= regenerate.(G.loop)
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
