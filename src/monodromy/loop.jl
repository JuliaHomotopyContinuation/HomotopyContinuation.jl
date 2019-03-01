export Triangle, Petal

#####################
# Nodes of the Loop
#####################

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

"""
    add!(node::Node, x; kwargs...)

Calls [`add!`](@ref) on the points of the Node.
"""
add!(node::Node, x; kwargs...) = add!(node.points, x; kwargs...)

"""
    iscontained(node::Node, x; kwargs...)

Calls [`iscontained`](@ref) on the points of the Node.
"""
function iscontained(node::Node, x; kwargs...)
    if node.points === nothing
        false
    else
        iscontained(node.points, x; kwargs...)
    end
end

"""
    unsafe_add!(node::Node, x)

Calls [`unsafe_add!`](@ref) on the points of the Node.
"""
function unsafe_add!(node::Node, x)
    if node.points !== nothing
        unsafe_add!(node.points, x)
    end
    nothing
end


#####################
# Edges of the Loop
#####################
"""
    Edge(start::Int, target::Int; usegamma=false)

Create an edge between two nodes references by the index. `usegamma` refers
to an optional weight of two random complex numbers γ=(γ₁, γ₀). These are the same
gammas as in [`ParameterHomotopy`](@ref).

    Edge(edge::Edge)

Construct an indentical `Edge` as `edge`, but with newly samples γ (if applicable).
"""
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
function Edge(edge::Edge)
    if edge.γ === nothing
        edge
    else
        γ = (randn(ComplexF64), randn(ComplexF64))
        Edge(edge.start, edge.target, γ)
    end
end

#######################
# Loop data structure
#######################
"""
    Loop(p::SVector, x::AbstractVector{<:Number}, nnodes::Int, options::MonodromyOptions; usegamma=true)

Construct a loop using `nnodes` nodes with parameters of the type of `p` and solutions
of the type of `x`. `usegamma` refers to the use weights on the edges. See also
[`Edge`](@ref).

    Loop(p::SVector, x::AbstractVector, nnodes::Int, options::MonodromyOptions; usegamma=true)

Construct the loop and add all points in `x` to it.

    Loop(style::LoopStyle, p::SVector, x::AbstractVector, options::MonodromyOptions)

Construct a loop using the defined style `style` with parameters of the type of `p` and solutions
of the type of `x`.
"""
struct Loop{N<:Node}
    nodes::Vector{N}
    edges::Vector{Edge}
end

function Loop(p₁::SVector, x₁::AbstractVector{<:Number}, nnodes::Int, options::MonodromyOptions; usegamma=true)
    n₁ = Node(p₁, x₁, is_main_node = true)
    nodes = [n₁]
    store_points = options.group_action_on_all_nodes && has_group_actions(options)
    for i = 2:nnodes
        p = options.parameter_sampler(p₁)
        push!(nodes, Node(p, x₁, store_points=store_points))
    end

    loop = map(i -> Edge(i - 1, i; usegamma=usegamma), 2:nnodes)
    push!(loop, Edge(nnodes, 1; usegamma=usegamma))

    Loop(nodes, loop)
end
function Loop(p₁::SVector, x::AbstractVector, nnodes::Int, options::MonodromyOptions; kwargs...)
    loop = Loop(p₁, first(x), nnodes, options; kwargs...)
    for xᵢ ∈ x
        add!(loop.nodes[1], xᵢ; tol=options.tol)
    end
    loop
end

"""
    solutions(loop::Loop)

Get the solutions of the loop.
"""
solutions(loop::Loop) = loop.nodes[1].points

"""
    nsolutions(loop::Loop)

Get the number solutions of the loop.
"""
nsolutions(loop::Loop) = length(solutions(loop))

"""
    mainnode(loop::Loop)

Get the main [`Node`](@ref) of the loop, i.e., the one with the original provided parameter.
"""
mainnode(loop::Loop) = loop.nodes[1]

nextedge(loop::Loop, edge::Edge) = loop.edges[edge.target]

##############
# Loop styles
##############
"""
    LoopStyle

Abstract type defining a style of a loop.
"""
abstract type LoopStyle end


"""
    Triangle(;weights=true)

A triangle is a loop consisting of the main node and two addtional nodes.
If `weights` is true the edges are equipped with additional random weights.
Note that this is usually only necessary for real parameters.
"""
struct Triangle <: LoopStyle
    useweights::Bool
end
Triangle(;useweights=true) = Triangle(useweights)

function Loop(strategy::Triangle, p::SVector, x::AbstractVector, options::MonodromyOptions)
    Loop(p, x, 3, options, usegamma=strategy.useweights)
end

"""
    Petal()

A petal is a loop consisting of the main node and one other node connected
by two edges with different random weights.
"""
struct Petal <: LoopStyle end
function Loop(strategy::Petal, p::SVector, x::AbstractVector, options::MonodromyOptions)
    Loop(p, x, 2, options, usegamma=true)
end


"""
    regenerate!(loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics)

Regenerate all random parameters in the loop in order to introduce a new monodromy action.
"""
function regenerate!(loop::Loop, options::MonodromyOptions, stats::MonodromyStatistics)
    main = mainnode(loop)

    # The first node is the main node and doesn't get touched
    for i ∈ 2:length(loop.nodes)
        loop.nodes[i] = Node(options.parameter_sampler(main.p), loop.nodes[i])
    end
    loop.edges .= Edge.(loop.edges)
    generated_parameters!(stats, length(main.points)) # bookkeeping
end


"""
    track(tracker, x::AbstractVector, edge::Edge, loop::Loop, stats::MonodromyStatistics)

Track `x` along the edge `edge` in the loop `loop` using `tracker`. Record statistics
in `stats`.
"""
function track(tracker::PathTracker, x::AbstractVector, edge::Edge, loop::Loop, stats::MonodromyStatistics)
    set_parameters!(tracker, edge, loop)
    track(tracker, x, stats)
end
function track(tracker::PathTracker, x::AbstractVector, stats::MonodromyStatistics)
    retcode = track!(tracker, x, 1.0, 0.0)
    trackedpath!(stats, retcode)
    retcode
end

"""x
    set_parameters!(tracker::PathTracker, e::Edge, loop::Loop)

Setup the parameters in the ParameterHomotopy in `tracker` to fit the edge `e`.
"""
function set_parameters!(tracker::PathTracker, e::Edge, loop::Loop)
    H = basehomotopy(tracker.homotopy)
    if !(H isa ParameterHomotopy)
        error("Base homotopy is not a ParameterHomotopy")
    end
    p = (loop.nodes[e.start].p, loop.nodes[e.target].p)
    set_parameters!(H, p, e.γ)
end
