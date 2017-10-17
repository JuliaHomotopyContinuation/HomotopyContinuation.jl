struct PathResult{T}
    retcode::Symbol
    solution::T

    residual::Float64
    newton_residual::Float64
    condition_jacobian::Float64
    windingnumber::Int
    homogenous_coordinate_magnitude::Float64

    startvalue::T

    iterations::Int
    endgame_iterations::Int
    npredictions::Int
end


struct Result{T} <: AbstractVector{PathResult{T}}
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
