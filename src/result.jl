export Result

struct PathResult{T}
    returncode::Symbol
    solution::Vector{T}
    singular::Bool

    residual::Float64
    newton_residual::Float64
    log10_condition_number::Float64
    windingnumber::Int
    angle_to_infinity::Float64
    real_solution::Bool

    startvalue::Vector{T}

    iterations::Int
    endgame_iterations::Int
    npredictions::Int
end

function PathResult(startvalue::AbstractVector, trackedpath_result::PathtrackerResult{T}, endgamer_result::EndgamerResult, solver::Solver) where T
    @unpack pathtracker, options = solver
    @unpack at_infinity_tol, singular_tol, abstol = options

    @unpack returncode, solution, windingnumber = endgamer_result

    if returncode == :success
        returncode = :isolated
    end

    residual, newton_residual, condition_number = residual_estimates(solution, pathtracker)

    # check whether startvalue was affine and our solution is projective
    N = length(startvalue)
    if length(solution) == N + 1
        # make affine

        homog_var = solution[1]
        affine_var = solution[2:end]
        angle_to_infinity = atan(abs(homog_var)/norm(affine_var))
        if angle_to_infinity < at_infinity_tol
            returncode = :at_infinity
            a, index = findmax(abs2.(affine_var))
            scale!(solution, inv(affine_var[index]))
        else
            scale!(affine_var, inv(homog_var))
            solution = affine_var
        end
    else
        angle_to_infinity = NaN
    end

    singular = windingnumber > 1 || condition_number > singular_tol

    if norm(imag(solution)) < condition_number * abstol
        real_solution = true
    else
        real_solution = false
    end

    PathResult{T}(
        returncode,
        solution,
        singular,
        residual,
        newton_residual,
        log10(condition_number),
        windingnumber,
        angle_to_infinity,
        real_solution,
        convert(typeof(solution), startvalue),
        trackedpath_result.iterations,
        endgamer_result.iterations,
        endgamer_result.npredictions
        )
end


function residual_estimates(solution, tracker::Pathtracker{Low}) where Low
    @unpack H, cfg, cache = tracker.low
    residual, newton_residual = residuals(H, solution, 0.0, cfg, cache)

    condition_number::Float64 = Homotopy.κ(H, solution, 0.0, cfg)

    residual, newton_residual, condition_number
end


function Base.show(io::IO, r::PathResult)
    println(io, typeof(r), ":")
    println(io, "* returncode: $(r.returncode)")
    println(io, "* solution: $(r.solution)")
    println(io, "* singular: $(r.singular)")
    println(io, "---------------------------------------------")
    println(io, "* iterations: $(r.iterations)")
    println(io, "* endgame iterations: $(r.endgame_iterations)")
    println(io, "* npredictions: $(r.npredictions)")
    println(io, "---------------------------------------------")
    println(io, "* newton_residual: $(@sprintf "%.3e" r.newton_residual)")
    println(io, "* residual: $(@sprintf "%.3e" r.residual)")
    println(io, "* log10 of the condition_number: $(@sprintf "%.3e" r.log10_condition_number)")
    println(io, "* windingnumber: $(r.windingnumber)")
    println(io, "* angle to infinity: $(round(r.angle_to_infinity,3))")
    println(io, "* real solution: $(r.real_solution)")
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


function Base.show(io::IO, result::Result{T}) where T
    println(io, typeof(result), ":")
    println(io, "* Total number of paths: $(length(result.pathresults))")
    println(io, "* Number of successfull paths: $(sum(r -> r.returncode == :isolated ? 1 : 0, result.pathresults))")
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
            Juno.render(i, Text("# isolated solutions → $(sum(r -> r.returncode == :isolated, s.pathresults))")),
            Juno.render(i, Text("# solutions at infinity → $(sum(r -> r.returncode == :at_infinity, s.pathresults))")),
            pathresults_t]
        return t
    end

end
