export Result

struct PathResult{T}
    returncode::Symbol
    solution::Vector{T}

    residual::Float64
    newton_residual::Float64
    condition_jacobian::Float64
    windingnumber::Int
    angle_to_infinity::Float64

    startvalue::Vector{T}

    iterations::Int
    endgame_iterations::Int
    npredictions::Int
end

function Base.show(io::IO, r::PathResult)
    println(io, typeof(r), ":")
    println(io, "* retcode: $(r.retcode)")
    println(io, "* solution: $(r.solution)")
    println(io, "---------------------------------------------")
    println(io, "* iterations: $(r.iterations)")
    println(io, "* endgame_iterations: $(r.endgame_iterations)")
    println(io, "* npredictions: $(r.npredictions)")
    println(io, "---------------------------------------------")
    println(io, "* newton_residual: $(r.newton_residual)")
    println(io, "* residual: $(r.residual)")
    println(io, "* condition_jacobian: $(r.condition_jacobian)")
    println(io, "* windingnumber: $(r.windingnumber)")
    println(io, "* homogenous_coordinate_magnitude: $(r.homogenous_coordinate_magnitude)")
end


mutable struct Result{T} <: AbstractVector{PathResult{T}}
    pathresults::Vector{PathResult{T}}
end

Base.start(result::Result) = start(result.pathresults)
Base.next(result::Result, state) = next(result.pathresults, state)
Base.done(result::Result, state) = done(result.pathresults, state)
Base.eltype(result::Result) = eltype(result.pathresults)
Base.length(result::Result) = length(result.pathresults)
Base.getindex(result::Result, i) = getindex(result.pathresults, i)
Base.setindex!(result::Result, x, i) = setindex!(result.pathresults, x, i)
Base.endof(result::Result) = endof(result.pathresults)
Base.size(result::Result, dim) = size(result.pathresults, dim)
Base.size(result::Result) = size(result.pathresults)


function Base.show(io::IO, result::Result{T}) where T
    println(io, typeof(result), ":")
    println(io, "* Total number of paths: $(length(result.pathresults))")
    println(io, "* Number of successfull paths: $(sum(r -> r.retcode == :success ? 1 : 0, result.pathresults))")
end


@require Juno begin
    using Juno
    function Juno.render(i::Juno.Inline, s::PathResult)
        t = Juno.render(i, Juno.defaultrepr(s))
        return t
    end

    function Juno.render(i::Juno.Inline, s::Result)
        t = Juno.render(i, Juno.defaultrepr(s))
        pathresults_t = Juno.render(i, s.pathresults)
        t[:children] = [
            Juno.render(i, Text("Total number of paths → $(length(s.pathresults))")),
            Juno.render(i, Text("Number of successfull paths → $(sum(r -> r.retcode == :success ? 1 : 0, s.pathresults))")),
            pathresults_t]
        return t
    end

end
