export endgame!


function endgame!(endgamer::Endgamer, x::AbstractVector, t)
    reset!(endgamer, x, t)
    endgame!(endgamer)
    endgamer
end

function endgame!(endgamer::Endgamer)
    # TODO: Is the iterator fast enough to just is it directly?
    start(endgamer)
    is_done = done(endgamer, 0)
    while !is_done
        next(endgamer, 0)
        is_done = done(endgamer, 0)
    end
end


Base.eltype(::T) where {T<:Endgamer} = T

Base.start(::Endgamer) = 0
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

@inline function Base.next(endgamer::Endgamer, state)
    predict!(endgamer, endgamer.cache)
    check_convergence!(endgamer)
    if endgamer.status != Failed
        moveforward!(endgamer)
    end
    endgamer, state + 1
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

    prednorm = projectivenorm2(p, pprev)
    if prednorm < options.abstol * options.abstol
        endgamer.status = Successfull
    elseif npredictions > 2
        pprev2 = predictions[end - 2]
        # if we do not half our error we are satisfied with a less strict tolerance
        if projectivenorm2(pprev, pprev2) / prednorm < 2
            if prednorm < options.abstol
                endgamer.status = Successfull
            end
        end
    end
end


@inline function moveforward!(endgamer::Endgamer)
    @unpack R, tracker, xs = endgamer
    λ = endgamer.options.geometric_series_factor

    endgamer.iter += 1

    λR = λ * R
    # we should track with high precision if necessary
    track!(tracker, last(xs), R, λR)
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
