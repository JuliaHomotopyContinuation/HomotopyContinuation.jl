export Result, gammatrick_gamma, PathResult

"""
    PathResult(startvalue, pathtracker_result, endgamer_result, solver)

Construct a `PathResult` for a given `startvalue`. `pathtracker_result` is the
[`PathtrackerResult`](@ref) until the endgame radius is reached. `endgamer_result`
is the [`EndgamerResult`](@ref) resulting from the corresponding endgame.

A `PathResult` contains:
* `returncode`: One of `:success`, `:at_infinity` or any error code from the `EndgamerResult`
* `solution::Vector{T}`: The solution vector. If the algorithm computed in projective space
and the solution is at infinity then the projective solution is given. Otherwise
an affine solution is given if the startvalue was affine and a projective solution
is given if the startvalue was projective.
* `residual::Float64`: The value of the infinity norm of `H(solution, 0)`.
* `newton_residual`: The value of the 2-norm of ``J_H(\\text{solution})^{-1}H(\\text{solution}, 0)``
* `log10_condition_number`: A high condition number indicates singularty. See [`Homotopy.κ`](@ref) for details.
    The value is the logarithmic condition number (with base 10).
* `windingnumber`: The estimated winding number
* `angle_to_infinity`: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
* `real_solution`: Indicates whether the solution is real given the defined tolerance `at_infinity_tol` (from the solver options).
* `startvalue`: The startvalue of the path
* `iterations`: The number of iterations the pathtracker needed.
* `endgame_iterations`: The number of steps in the geometric series the endgamer did.
* `npredictions`: The number of predictions the endgamer did.
* `predictions`: The predictions of the endgamer.
"""
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
    predictions::Vector{Vector{T}}
end

function PathResult(startvalue::AbstractVector, trackedpath_result::PathtrackerResult{T}, endgamer_result::EndgamerResult, solver::Solver) where T
    @unpack pathtracker, options = solver
    @unpack at_infinity_tol, singular_tol, abstol = options

    @unpack returncode, solution, windingnumber = endgamer_result


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
        endgamer_result.npredictions,
        endgamer_result.predictions
        )
end


function residual_estimates(x, tracker::Pathtracker{Low}) where Low
    @unpack H, cfg, cache = tracker.low

    res = evaluate(H, x, 0.0, cfg)
    jacobian = Homotopy.jacobian(H, x, 0.0, cfg, true)
    residual = norm(res, Inf)
    newton_residual::Float64 = norm(pinv(jacobian) * res)

    condition_number::Float64 = Homotopy.κ(H, x, 0.0, cfg)

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

"""
    Result(pathresults, solver)

A thin wrapper around the [`PathResult`](@ref)s of the `Solver` instance.
`Result` behaves like an array of `PathResult`s but also contains some additional
informations. For example you can obtain the γ which was used for the gammatrick.
"""
struct Result{T} <: AbstractVector{PathResult{T}}
    gamma::Complex128
    pathresults::Vector{PathResult{T}}
end

Result(pathresults, solver) = Result(solver.gamma, pathresults)

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

"""
    gammatrick_gamma(result)

Return the γ which was used for the gammatrick.
"""
gammatrick_gamma(r::Result) = r.gamma

function Base.show(io::IO, result::Result{T}) where T
    println(io, typeof(result), ":")
    println(io, "# paths: $(length(result.pathresults))")
    println(io, "# successfull paths: $(sum(r -> r.returncode == :success ? 1 : 0, result.pathresults))")
    println(io, "# solutions at infinity → $(sum(r -> r.returncode == :at_infinity, s.pathresults))")
    println(io, "# singular solutions → $(sum(r -> r.singular, s.pathresults))")
    println(io, "# real solutions → $(sum(r -> r.real_solution, s.pathresults))")
    # println(io, "gamma → $(result.gamma)")
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
            Juno.render(i, Text("# paths → $(length(s.pathresults))")),
            Juno.render(i, Text("# successfull paths → $(sum(r -> r.returncode == :success, s.pathresults))")),
            Juno.render(i, Text("# solutions at infinity → $(sum(r -> r.returncode == :at_infinity, s.pathresults))")),
            Juno.render(i, Text("# singular solutions → $(sum(r -> r.singular, s.pathresults))")),
            Juno.render(i, Text("# real solutions → $(sum(r -> r.real_solution, s.pathresults))")),
            # Juno.render(i, Text("gamma → $(s.gamma)")),
            pathresults_t]
        return t
    end

end
