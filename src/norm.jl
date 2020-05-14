export AbstractNorm, InfNorm, WeightedNorm, distance, norm, init!, update!

using Base: @propagate_inbounds
import LinearAlgebra: norm

"""
    AbstractNorm

An `AbstractNorm` represents any norm of a vector space.
All norms are callable. `norm(x)` computes the norm of `x` and `norm(x,y)` computes the distance
`norm(x - y).`
"""
abstract type AbstractNorm end

"""
    distance(u, v, norm::AbstractNorm)

Compute the distance ||u-v|| with respect to the given norm `norm`.
"""
distance(u, v, norm::AbstractNorm) = MethodError(distance, (u, v, norm))

"""
    norm(u, norm::AbstractNorm)

Compute the norm ||u|| with respect to the given norm `norm`.
"""
LinearAlgebra.norm(u, norm::AbstractNorm) = MethodError(LinearAlgebra.norm, (u, norm))

##################
## WeightedNorm ##
##################

Base.@kwdef struct WeightedNormOptions
    scale_min::Float64 = sqrt(eps())
    scale_abs_min::Float64 = min(scale_min^2, sqrt(eps()))
    scale_max::Float64 = exp2(div(1023, 2))
end

Base.show(io::IO, opts::WeightedNormOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::WeightedNormOptions) = opts

"""
    WeightedNorm(d::AbstractVector, norm::AbstractNorm; options...)
    WeightedNorm(norm::AbstractNorm, n::Integer; options...)
    WeightedNorm(norm::AbstractNorm, x::AbstractVector; options...)

A `WeightedNorm` represents a weighted variant of norm `norm` of a `n`-dimensional vector space.`
A norm ``||x||`` is weighted by introducing a vector of additional weights `d` such that
the new norm is ``||D⁻¹x||`` where ``D`` is the diagonal matrix with diagonal ``d`.
The `WeightedNorm` is desigened to change the weights dynamically by using [`init!(::WeightedNorm, x)`](@ref)
and [`update!(::WeightedNorm, x)`](@ref). The weights are there constructed such that
``||D⁻¹x|| ≈ 1.0``.
The weights can be accessed and changed by indexing.

### Options

* `scale_min = sqrt(eps())`: The minimal size of `dᵢ` is `scale_min` time the (weighted) norm of `x`.
* `scale_abs_min = min(scale_min^2, 200 * sqrt(eps()))`: The absolute minimal size of `dᵢ`.
* `scale_max = 1.0 / eps() / sqrt(2)`: The absolute maximal size of `dᵢ`

"""
struct WeightedNorm{N<:AbstractNorm} <: AbstractNorm
    weights::Vector{Float64}
    norm::N
    options::WeightedNormOptions
end
function WeightedNorm(
    weights::Vector{Float64},
    norm::AbstractNorm;
    scale_min = sqrt(eps()),
    scale_abs_min = min(scale_min^2, sqrt(eps())),
    scale_max = exp2(div(1023, 2)),
)
    opts = WeightedNormOptions(scale_min, scale_abs_min, scale_max)
    WeightedNorm(weights, norm, opts)
end
WeightedNorm(norm::AbstractNorm, x::AbstractVector; opts...) =
    WeightedNorm(norm, length(x); opts...)
WeightedNorm(norm::AbstractNorm, n::Integer; opts...) = WeightedNorm(ones(n), norm; opts...)

(N::WeightedNorm)(x::AbstractVector) = norm(x, N)
(N::WeightedNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    weights(WN::WeightedNorm)

Returns the weights of the weighted norm.
"""
weights(WN::WeightedNorm) = WN.weights

Base.length(WN::WeightedNorm) = length(WN.weights)
Base.size(WN::WeightedNorm) = (length(WN),)
Base.copyto!(WN::WeightedNorm, x) = copyto!(WN.weights, x)
@propagate_inbounds Base.getindex(WN::WeightedNorm, i::Integer) = getindex(WN.weights, i)
@propagate_inbounds Base.setindex!(WN::WeightedNorm, x, i::Integer) =
    setindex!(WN.weights, x, i)


"""
    init!(w::WeightedNorm, x::AbstractVector)

Setup the weighted norm `w` for `x`.
"""
function init!(w::WeightedNorm, x::AbstractVector)
    point_norm = w.norm(x)
    for i = 1:length(w)
        wᵢ = fast_abs(x[i])
        if wᵢ < w.options.scale_min * point_norm
            wᵢ = w.options.scale_min * point_norm
        elseif wᵢ > w.options.scale_max * point_norm
            wᵢ = w.options.scale_max * point_norm
        end
        w[i] = max(wᵢ, w.options.scale_abs_min)
    end
    w
end
init!(n::AbstractNorm, ::AbstractVector) = n

"""
    update!(w::WeightedNorm, x::AbstractVector)

Update the weighted norm `w` for `x`, this will interpolate between the previous weights
and the norm of `x`.
"""
function update!(w::WeightedNorm, x::AbstractVector)
    norm_x = w(x)
    for i = 1:length(x)
        # wᵢ = fast_abs(x[i])#
        wᵢ = (fast_abs(x[i]) + w[i]) / 2
        if wᵢ < w.options.scale_min * norm_x
            wᵢ = w.options.scale_min * norm_x
        elseif wᵢ > w.options.scale_max * norm_x
            wᵢ = w.options.scale_max * norm_x
        end
        if isfinite(wᵢ)
            w[i] = max(wᵢ, w.options.scale_abs_min)
        end
    end
    w
end
update!(n::AbstractNorm, ::AbstractVector) = n


####################
## Specific norms ##
####################
"""
    InfNorm <: AbstractNorm

The infinity or [maximum norm](https://en.wikipedia.org/wiki/Norm_(mathematics)#Maximum_norm).
"""
struct InfNorm <: AbstractNorm end

function distance(x::AbstractVector, y::AbstractVector, ::InfNorm)
    n = length(x)
    @boundscheck n == length(y)
    @inbounds dmax = abs2(x[1] - y[1])
    for i = 2:n
        @inbounds dᵢ = abs2(x[i] - y[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    d = sqrt(dmax)
    if isinf(d)
        @inbounds dmax = abs(x[1] - y[1])
        for i = 2:n
            @inbounds dᵢ = abs2(x[i] - y[i])
            dmax = max(dmax, dᵢ)
        end
        return dmax
    else
        return d
    end
end
function distance(x::AbstractVector, y::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2((x[1] - y[1]) / w[1])
    for i = 2:n
        @inbounds dᵢ = abs2((x[i] - y[i]) / w[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    d = sqrt(dmax)
    if isinf(d)
        @inbounds dmax = abs((x[1] - y[1]) / w[1])
        for i = 2:n
            @inbounds dᵢ = abs((x[i] - y[i]) / w[i])
            dmax = max(dmax, dᵢ)
        end
        return dmax
    else
        return d
    end
end

function norm(x::AbstractVector, ::InfNorm)
    n = length(x)
    @inbounds dmax = abs2(x[1])
    for i = 2:n
        @inbounds dᵢ = abs2(x[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    d = sqrt(dmax)

    if isinf(d)
        @inbounds dmax = abs(x[1])
        for i = 2:n
            @inbounds dᵢ = abs(x[i])
            dmax = max(dmax, dᵢ)
        end
        return dmax
    else
        return d
    end
end
function norm(x::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2(x[1] / w[1])
    for i = 2:n
        @inbounds dᵢ = abs2(x[i] / w[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    d = sqrt(dmax)

    if isinf(d)
        @inbounds dmax = abs(x[1] / w[1])
        for i = 2:n
            @inbounds dᵢ = abs(x[i] / w[i])
            dmax = max(dmax, dᵢ)
        end
        return dmax
    else
        return d
    end
end
(N::InfNorm)(x::AbstractVector) = norm(x, N)
(N::InfNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)
