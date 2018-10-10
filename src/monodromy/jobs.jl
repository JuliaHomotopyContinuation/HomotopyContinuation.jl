struct Job{N, T}
    x::SVector{N, T}
    edge::Edge
end

struct JobResult{N, T}
    retcode::Symbol
    x::SVector{N, T}
end


function execute(job::Job, G::Graph, tracker, stats::Statistics)
    y, retcode = track(tracker, job.x, job.edge, G, stats)
    JobResult(retcode, y)
end

function process!(queue::Vector{<:Job}, job::Job, res::JobResult, G::Graph, options::Options)
    if res.retcode ≠ :success
        return :incomplete
    end
    node = nodes(G)[job.edge.target]
    if add!(node, res.x, tol=options.tol)
        if isdone(node, res.x, options)
            return :done
        end
        next_edge = nextedge(G, job.edge)
        push!(queue, Job(res.x, next_edge))

        # handle group actions
        if node.apply_group_action
            for yᵢ in options.group_actions(res.x)
                if add!(node, yᵢ, tol=options.tol)
                    if isdone(node, yᵢ, options)
                        return :done
                    end
                    push!(queue, Job(yᵢ, next_edge))
                end
            end
        end
    end
    return :incomplete
end

function isdone(node::Node, x, options::Options)
    !node.main_node && return false

    options.done_callback(x) ||
    length(node.points) ≥ options.target_solutions_count
end
