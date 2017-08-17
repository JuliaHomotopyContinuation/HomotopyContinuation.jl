"""
    PathResult(result, returncode, startvalue, iterations, steps, trace)

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

function Base.show(io::IO, res::PathResult)
    println(io, typeof(res),":")
    println(io, "* result: ", result)
    println(io, "* returncode: ", returncode)
    println(io, "* startvalue: ", startvalue)
    println(io, "* iterations: ", iterations)
    println(io, "* laststep: ", laststep)
    println(io, "* steps: ", length(steps), " entries")
    println(io, "* trace: ", length(trace), " entries")
end



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
    algorithm::AbstractPredictorCorrectorHomConAlgorithm,
    start::S,
    finish::S;
    maxiterations=10000,
    tolerance=1e-4,
    refinement_tolerance=1e-8,
    initial_steplength=0.01,
    successfull_steps_until_steplength_increase=3,
    steplength_increase_factor=2.0,
    steplength_decrease_factor=0.5,
    record_trace=false,
    record_steps=false,
    correction_step_maxiterations=3,
    debug=false
) where {S<:Number,T<:Number}

    # we only need to compute these once
    J_H = differentiate(H)
    ∂H∂t = ∂t(H)

    # We want to track a path z(t) from z(start) to z(finish).
    # We reformulate it as a problem z(γ(s)) where γ: [1, 0] -> line from start to finish
    γlength = finish - start
    γ(s) = start + (1 - s) * γlength

    # some setup before we start the main loop
    x = copy(startvalue)
    u = similar(x)

    steps::Vector{S} = [];
    trace::Vector{Vector{T}} = [];

    s = 1.0
    steplength = initial_steplength
    sucessive_successes = 0
    k = 0

    while s > 0 && k < maxiterations
        if record_steps
            push!(steps, s)
        end
        if record_trace
            push!(trace, x)
        end

        Δs= min(steplength, s)

        u .= predict(algorithm, H, J_H, ∂H∂t, x, γ(s), Δs * γlength)

        if debug
            println("")
            println("t: $t Δt: $Δt, iteration: $k")
            println("x: $x, norm: $(norm(evaluate(H,x,t)))")
            println("predicted: $x, norm: $(norm(evaluate(H,x,t)))")
        end

        converged = correct!(u, algorithm, H, J_H, x, γ(s - Δs), tolerance, correction_step_maxiterations)

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
        refinement_converged = correct!(u, algorithm, H, J_H, x, γ(0.0), refinement_tolerance, 10)
        if refinement_converged
            x .= u
        end
    end

    if k >= maxiterations
        retcode = :MaxIterations
    elseif norm(evaluate(H, x, γ(0.0))) > tolerance
        retcode = :Diverged
    else
        retcode = :Success
    end

    PathResult(x, retcode, startvalue, k, γ(s), steps, trace)
end
