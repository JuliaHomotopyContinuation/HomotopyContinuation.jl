export solve

"""
    solve(H::AbstractHomotopy{T}, startvalue::Vector{T}, [algorithm,] endgame_start=0.1, kwargs...)

Track the path ``x(t)`` implicitly defined by ``H(x(t),t)`` from `t=1` to `t=0` where
`x(1)=startvalue`. Returns a [`Result`](@ref) instance.

It uses the defined `algorihm` for the prediction and correction step.
If no `algorithm` is defined the in general best performing algorithm is choosen. Currently
this is `Spherical`.

## Optional arguments
* `endgame_start=0.1`: When to switch to the Cauchy endgame, for a value of t=0.0 no endgame happens.
*  optional arguments to [`pathtracking`](@ref)
*  optional arguments to [`cauchyendgame`](@ref)

    solve(H::AbstractHomotopy, startvalues, algorithm, report_progess=false, kwargs...)

Track multiple paths with the given `startvalues`. Returns a vector of [`Result`](@ref) instances,
one per start value.
If `report_progess=true` the current progess will be logged.
"""
#
# Most general loop for multiple start values
#
function solve(
    H::AbstractHomotopy{T},
    startvalues,
    algorithm::APCA;
    report_progress=false,
    kwargs...
) where {T<:Number}
    if eltype(startvalues) != Vector{T}
        error("Expected as start values an iterable with elements of type Vector{$T} " *
            "but got $(eltype(startvalues))")
    end

    H = homogenize_if_necessary(H, algorithm)

    J_H! = Homotopy.jacobian!(H)
    Hdt! = Homotopy.dt!(H)

    if report_progress
        println("Total number of paths to track: $length(startvalues)")
    end

    map(1:length(startvalues), startvalues) do index, startvalue
        if report_progress
            println("Start to track path $index")
        end

        result = solve(H, J_H!, Hdt!, startvalue, algorithm; kwargs...)

        if report_progress
            println("Path $index returned: $(result.returncode)")
        end

        result
    end
end

# Here are now the actual solve routines for a single path
# (different implementations for affine and projective
# since there is currently no endgame for affine algorithms)
function solve(
    H::AbstractHomotopy{T},
    J_H!::Function,
    Hdt!::Function,
    startvalue::Vector{T},
    algorithm::APCA{Val{false}};
    start::Float64=1.0,
    trackpathkwargs...
) where {T<:Number}
    @assert length(startvalue) == Homotopy.nvariables(H)
        "A start_value has length $(length(startvalue)). Excepted length $(length(Homotopy.nvariables(H)))."

    pathresult = trackpath(H, J_H!, Hdt!, startvalue, algorithm, start, 0.0; trackpathkwargs...)
    Result(pathresult, startvalue, algorithm)
end

function solve(
    H::AbstractHomotopy{T},
    J_H!::Function,
    Hdt!::Function,
    startvalue::Vector{T},
    algorithm::APCA{Val{true}};
    start::Float64=1.0,
    endgame_start::Float64=0.1,
    tolerance_infinity=1e-6,
    kwargs...
) where {T<:Number}
    @assert Homotopy.nvariables(H) == length(H) + 1 "Expected a projective homotopy, got $(H)"

    x = embed_projective_if_necessary(startvalue, H)
    # we have to split kwargs
    tpkwargs = trackpathkwargs(kwargs)
    ckwargs = cauchykwargs(kwargs)

    pathresult = trackpath(H, J_H!, Hdt!, x, algorithm, start, endgame_start; tpkwargs...)

    if !issuccessfull(pathresult) || (endgame_start <= 0.0)
        return Result(pathresult, startvalue, algorithm, tolerance_infinity)
    end

    endgameresult =
        cauchyendgame(
            H, J_H!, Hdt!, result(pathresult), endgame_start, algorithm;
            tolerance_infinity=tolerance_infinity,
            ckwargs...,
            tpkwargs...
        )
    Result(pathresult, endgameresult, startvalue, algorithm, tolerance_infinity)
end

homogenize_if_necessary(H::AbstractHomotopy, alg::APCA{Val{true}}) = Homotopy.homogenize(H)
homogenize_if_necessary(H::AbstractHomotopy, alg::APCA{Val{false}}) = H
"""
    embed_projective_if_necessary(x, H, ::AbstractPredictorCorrectorAlgorithm{true})

Embeds a vector into the projective space if necessary, i.e. if it's length is one less
than the number of variables of `H`. `H` is assumed to be homogenized. After the (eventual)
embedding the value is normalized.
"""
function embed_projective_if_necessary(x::Vector{T}, H::AbstractHomotopy{T}) where T
    N = Homotopy.nvariables(H)
    n = length(x)
    if N == n
        return x
    elseif N - 1 == n
        return [one(T); x]
    else
        return error("A start value has length $n. Excepted length $N or $(N-1).")
    end
end
