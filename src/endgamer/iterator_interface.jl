Base.start(::Endgamer) = 0

@inline function Base.next(endgamer::Endgamer, state)
    predict!(endgamer, endgamer.cache)
    check_convergence!(endgamer)
    if endgamer.status != Failed
        moveforward!(endgamer)
    end
    endgamer, state + 1
end

@inline function Base.done(endgamer::Endgamer, state)
    @unpack status = endgamer
    if status == Successfull || status == Failed
        return true
    end

    if endgamer.R ≤ eps(Float64)
        endgamer.status = Failed
        endgamer.failurecode = :ill_conditioned_zone
        return true
    end

    false
end


@inline function moveforward!(endgamer::Endgamer)
    @unpack R, tracker, xs = endgamer
    λ = endgamer.options.geometric_series_factor

    endgamer.iter += 1

    λR = λ * R
    # we should track with high precision if necessary
    run!(tracker, last(xs), R, λR)
    retcode, sol = solution(tracker)

    if retcode != :success
        endgamer.status = Failed
        endgamer.failurecode = :ill_conditioned_zone
        return nothing
    end

    push!(endgamer.xs, sol)
    endgamer.R = λR
    nothing
end

@inline function check_convergence!(endgamer::Endgamer)
    @unpack predictions, options = endgamer

    npredictions = length(predictions)
    # we need at least 2 predictions
    if npredictions < 2
        return nothing
    end
    if npredictions > 1
        endgamer.status = Started
    end

    p, pprev = predictions[end], predictions[end-1]

    if sqrt(projectivenorm2(p, pprev)) < options.abstol
        endgamer.status = Successfull
    end
end


function run!(endgamer::Endgamer)
    # TODO: Is the iterator fast enough to just is it directly?
    while endgamer.status != Successfull && endgamer.status != Failed
        if endgamer.R ≤ eps(Float64)
            endgamer.status = Failed
            endgamer.failurecode = :ill_conditioned_zone
            break
        end
        predict!(endgamer, endgamer.cache)
        check_convergence!(endgamer)
        if endgamer.status != Failed
            moveforward!(endgamer)
        end
    end
    endgamer
end

function run!(endgamer::Endgamer, x, t)
    reset!(endgamer, x, t)
    run!(endgamer)
    endgamer
end
