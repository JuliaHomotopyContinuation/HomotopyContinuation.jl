"""
    PathResult(result, returncode, startvalue, iterations, laststep, steps, trace)

Construct a result of `trackpath`.

## Fields
* `result::Vector{T}`
* `returncode::Symbol`: :Success or :MaxIterations or :Diverged
* `startvalue::Vector{T}`
* `iterations::Int`
* `laststep::S`: If :Success this is just `finish`
* `steps::Vector{S}`: Empty if `record_steps=false` (default)
* `trace::Vector{Vector{T}}`: Empty if `record_trace=false` (default)
"""
struct PathResult{T<:Number, S<:Number}
    result::Vector{T}
    returncode::Symbol

    startvalue::Vector{T}
    iterations::Int
    laststep::S
    steps::Vector{S}
    trace::Vector{Vector{T}}
end

"""
    CauchyEndgameResult(result, returncode, iterations, R, λ, trace)

Construct a result of a `cauchyendgame`.

## Fields
* `result::Vector{T}`
* `returncode::Symbol`: `:Success`, `:IllConditionedZone` (the path tracking failed)
or `:MachineEpsilon`: (`t` got smaller than machine epsilon)
* `iterations::Int`: How many iterations of the power series was done
* `R`: The radius of the endgame zone
* `λ`: The factor of the geometric series ``λ^kR``
* `trace::Vector{Vector{T}}`: The points ``x(λ^kR)`` for `k=0,1,...,iterations`
"""
struct CauchyEndgameResult{T<:Number}
    result::Vector{T}
    returncode::Symbol
    iterations::Int
    R::Float64
    λ::Float64
    trace::Vector{Vector{T}}
end

"""
    Result(solution, returncode, iterations, affine_solution, projective_solution, startvalue, steps, trace)

Construct a result of `solve` with an `AbstractPredictorCorrectorAlgorithm`.

## Fields
* `solution::Vector{T}`: This has the same length as startvalue
* `returncode::Symbol`: :Success, :MaxIterations, :Diverged or :AtInfinity
* `iterations::Int`
* `endgame_iterations::Int`
* `affine_solution::Vector{T}`
* `project_solution::Vector{T}`
* `startvalue::Vector{T}`
* `trace::Vector{Vector{T}}`: Empty if `record_trace=false` (default)
"""
struct Result{T<:Number}
    solution::Vector{T}
    returncode::Symbol
    iterations::Int
    endgame_iterations::Int

    affine_solution::Vector{T}
    projective_solution::Vector{T}

    startvalue::Vector{T}
    #* `steps::Vector{S}`: Empty if `record_steps=false` (default)
    #steps::Vector{Float64}
    trace::Vector{Vector{T}}
end

function Base.show(io::IO, ::MIME"text/plain", res::Result)
    println(io, typeof(res),":")
    println(io, "------------------------------")
    println(io, "* solution: ", res.solution)
    println(io, "* returncode: ", res.returncode)
    println(io, "------------------------------")
    println(io, "* iterations: ", res.iterations)
    println(io, "* endgame_iterations: ", res.endgame_iterations)
    println(io, "* affine_solution: ", res.affine_solution)
    println(io, "* projective_solution: ", res.projective_solution)
    println(io, "* startvalue: ", res.startvalue)
    #println(io, "* steps: ", length(res.steps), " entries")
    println(io, "* trace: ", length(res.trace), " entries")
end

function Base.show(io::IO, ::MIME"text/plain", res::CauchyEndgameResult)
    println(io, typeof(res),":")
    println(io, "------------------------------")
    println(io, "* result: ", res.result)
    println(io, "* returncode: ", res.returncode)
    println(io, "------------------------------")
    println(io, "* iterations: ", res.iterations)
    println(io, "* R: ", res.R)
    println(io, "* λ: ", res.λ)
    println(io, "* trace: ", length(res.trace), " entries")
end

function Base.show(io::IO, ::MIME"text/plain", res::PathResult)
    println(io, typeof(res),":")
    println(io, "------------------------------")
    println(io, "* result: ", res.result)
    println(io, "* returncode: ", res.returncode)
    println(io, "------------------------------")
    println(io, "* startvalue: ", res.startvalue)
    println(io, "* iterations: ", res.iterations)
    println(io, "* laststep: ", res.laststep)
    println(io, "* steps: ", length(res.steps), " entries")
    println(io, "* trace: ", length(res.trace), " entries")
end

"""
    nosuccess(result)

Checks whether the returncode is something else than `:Success`
"""
nosuccess(result::PathResult) = result.returncode != :Success
nosuccess(result::Result) = result.returncode != :Success
nosuccess(result::CauchyEndgameResult) = result.returncode != :Success
