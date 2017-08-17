prepare_homotopy(H::AbstractHomotopy, alg::APCA{true}) = homogenize(H)
prepare_homotopy(H::AbstractHomotopy, alg::APCA{false}) = H

"""
    prepare_startvalue(H, startvalue, ::AbstractPredictorCorrectorAlgorithm{true})

Embeds a vector into the projective space if necessary, i.e. if it's length is one less
than the number of variables of `H`. `H` is assumed to be homogenized. After the (eventual)
embedding the value is normalized.
"""
function prepare_startvalue(H::AbstractHomotopy, startvalue, ::APCA{true})
    N = nvariables(H)
    n = length(startvalue)
    if N == n
        normalize(startvalue)
    elseif N - 1 == n
        normalize!([1; startvalue])
    else
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
end

function prepare_startvalue(H::AbstractHomotopy, startvalue, ::APCA{false})
    N = nvariables(H)
    n = length(startvalue)
    if N == n
        return startvalue
    else
        return error("A start_value has length $n. Excepted length $N.")
    end
end
