export trackpath, PathResult, result, returncode, iterations, startvalue, pathtrace, pathsteps, laststep, issuccessfull
"""
    trackpath(H::AbstractHomotopy, startvalue, algorithm::PredictorCorrector, start, finish; kwargs...)

Tracks the path ``z(t)`` which is implicitly defined via ``H(x,t)`` from
`z(start)=H(startvalue, start)` to `z(finish)= H(result, finish)` with the given
PredictorCorrector `algorithm`.

We reformulate the path ``z(t)`` as a path z(γ(s)) where γ: [1, 0] -> line from start to finish

## Optional arguments
* `maxiterations=10000`
* `tolerance=1e-4`: The tolerance used for the correction step
* `refinement_tolerance=1e-8`: The result will be refined to the given accuracy
* `initial_steplength=0.01`: The initial_steplength (w.r.t s)
* `correction_step_maxiterations=3`
* `successfull_steps_until_steplength_increase=3`
* `steplength_increase_factor=2.0`
* `steplength_decrease_factor=0.5`
* `record_trace=false`
* `record_steps=false`
* `debug=false`: Print *(quite a lot)* different values to debug behaviour.
"""
function trackpath(
    H::AbstractHomotopy{T},
    startvalue::Vector{T},
    algorithm::AbstractPredictorCorrectorAlgorithm,
    start::S,
    finish::S;
    kwargs...
) where {S<:Number,T<:Number}
    J_H! = Homotopy.jacobian!(H)
    Hdt! = Homotopy.dt!(H)
    trackpath(H, J_H!, Hdt!, startvalue, algorithm, start, finish; kwargs...)
end

function trackpath(
    H::AbstractHomotopy{T},
    J_H!::Function, # result of Homotopy.jacobian!, i.e. takes (U, x, t)
    Hdt!::Function, # result of Homotopy.dt!, i.e. takes (u, x, t)
    startvalue::Vector{T},
    algorithm::AbstractPredictorCorrectorAlgorithm,
    start::S,
    finish::S;
    # !keep kwargs in sync with `istrackpathkwarg` below!
    maxiterations=10000,
    tolerance=1e-4,
    refinement_tolerance=1e-12,
    initial_steplength=0.05,
    successfull_steps_until_steplength_increase=3,
    steplength_increase_factor=2.0,
    steplength_decrease_factor=0.5,
    record_trace=false,
    record_steps=false,
    correction_step_maxiterations=3,
    debug=false
) where {S<:Number,T<:Number}
    # We want to track a path z(t) from z(start) to z(finish).
    # We reformulate it as a problem z(γ(s)) where γ: [1, 0] -> line from start to finish
    γlength = finish - start
    γ(s) = start + (1 - s) * γlength

    # some setup before we start the main loop
    n = nvariables(H)
    x = copy(startvalue)
    u = similar(x)
    A = zeros(T, n, n)
    b = zeros(T, n)

    steps::Vector{S} = [];
    trace::Vector{Vector{T}} = [];

    s = 1.0
    steplength = initial_steplength
    sucessive_successes = 0
    k = 0

    while s > 0 && k < maxiterations
        if record_steps; push!(steps, s); end
        if record_trace; push!(trace, x); end

        Δs = min(steplength, s)
        predict!(u, A, b, H, J_H!, Hdt!, x, γ(s), Δs * γlength, algorithm)

        if debug
            println("")
            println("t: $t Δt: $Δt, iteration: $k")
            println("x: $x, ∞-norm: $(norm(evaluate(H,x,t), Inf))")
            println("predicted: $x, ∞-norm: $(norm(evaluate(H,x,t), Inf))")
        end

        converged = correct!(u, A, b, H, J_H!, x, γ(s - Δs), tolerance, correction_step_maxiterations, algorithm)

        if converged
            x .= u
            s -= Δs

            sucessive_successes += 1
            if sucessive_successes == successfull_steps_until_steplength_increase
                steplength *= steplength_increase_factor
                sucessive_successes = 0
            end
        else
            if debug
                println("correction failed!")
            end
            sucessive_successes = 0
            steplength *= steplength_decrease_factor
        end

        k += 1
    end

    if s ≈ 0 && norm(evaluate(H, x, γ(0.0))) < tolerance
        refinement_converged = correct!(u, A, b, H, J_H!, x, γ(0.0), refinement_tolerance, 10, algorithm)
        if refinement_converged
            x .= u
        end
    end

    if k >= maxiterations
        retcode = :MaxIterations
    elseif norm(evaluate(H, x, γ(0.0))) > refinement_tolerance
        retcode = :Diverged
    else
        retcode = :Success
    end

    PathResult(x, retcode, startvalue, k, γ(s), steps, trace)
end

function istrackpathkwarg(kwarg)
    symb = first(kwarg)
    symb == :maxiterations ||
    symb == :tolerance ||
    symb == :refinement_tolerance ||
    symb == :initial_steplength ||
    symb == :successfull_steps_until_steplength_increase ||
    symb == :steplength_increase_factor ||
    symb == :steplength_decrease_factor ||
    symb == :record_trace ||
    symb == :record_steps ||
    symb == :correction_step_maxiterations ||
    symb == :debug
end
trackpathkwargs(kwargs) = filter(istrackpathkwarg, kwargs)

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

result(r::PathResult) = r.result
returncode(r::PathResult) = r.returncode
iterations(r::PathResult) = r.iterations
startvalue(r::PathResult) = r.startvalue
laststep(r::PathResult) = r.laststep
pathtrace(r::PathResult) = r.trace
pathsteps(r::PathResult) = r.steps

issuccessfull(r::PathResult) = returncode(r) == :Success
