export euclidean_distance, euclidean_norm

"""
    AbstractNorm

An `AbstractNorm` represents any norm of a vector space.
"""
abstract type AbstractNorm end

"""
    distance(u, v, norm::AbstractNorm)

Compute the distance ||u-v|| with the given norm.
"""
function distance end


##################
## WeightedNorm ##
##################
"""
    WeightedNormOptions(;scale_min=0.01,
                     scale_abs_min=min(scale_min^2, 200 * sqrt(eps()),
                     scale_max=1.0 / eps() / sqrt(2)

Parameters for the auto scaling of the variables.
"""
struct WeightedNormOptions
    scale_min::Float64
    scale_abs_min::Float64
    scale_max::Float64
end
function WeightedNormOptions(;scale_min=0.001,
                             scale_abs_min=min(scale_min^2, 200 * sqrt(eps())),
                             scale_max=1.0 / eps() / sqrt(2))
    WeightedNormOptions(scale_min, scale_abs_min, scale_max)
end
Base.show(io::IO, opts::WeightedNormOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::WeightedNormOptions) = opts

"""
    WeightedNorm(norm, n::Integer)
    WeightedNorm(norm, x::AbstractVector)

A `WeightedNorm` represents a weighted variant of norm `norm` of a `n`-dimensional vector space.`
"""
struct WeightedNorm{N<:AbstractNorm} <: AbstractNorm
    weights::Vector{Float64}
    norm::N
    options::WeightedNormOptions
end
function WeightedNorm(weights::Vector{Float64}, norm::AbstractNorm; opts...)
    WeightedNorm(weights, norm, WeightedNormOptions(opts...))
end
WeightedNorm(norm::AbstractNorm, x::AbstractVector; opts...) = WeightedNorm(norm, length(x); opts...)
WeightedNorm(norm::AbstractNorm, n::Integer; opts...) = WeightedNorm(ones(n), norm; opts...)

(N::WeightedNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::WeightedNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    weights(WN::WeightedNorm)

Returns the inverse weights of the weighted norm.
"""
weights(WN::WeightedNorm) = WN.weights

Base.length(WN::WeightedNorm) = length(WN.weights)
Base.size(WN::WeightedNorm) = (length(WN),)
Base.copyto!(WN::WeightedNorm, x) = copyto!(WN.weights, x)
@propagate_inbounds Base.getindex(WN::WeightedNorm, i::Integer) = getindex(WN.weights, i)
@propagate_inbounds Base.setindex!(WN::WeightedNorm, x, i::Integer) = setindex!(WN.weights, x, i)


"""
    init!(w::WeightedNorm, x::AbstractVector)

Setup the weighted norm `w` for `x`.
"""
function init!(w::WeightedNorm, x::AbstractVector)
    point_norm = w.norm(x)
    for i in 1:length(w)
        wᵢ = abs(x[i])
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

Ipdate the weighted norm `w` for `x`, this will interpolate between the previous weights
and the norm of `x`.
"""
function update!(w::WeightedNorm, x::AbstractVector)
    norm_x = w(x)
    for i in 1:length(x)
        wᵢ = (abs(x[i]) + w[i]) / 2
        if wᵢ < w.options.scale_min * norm_x
            wᵢ = w.options.scale_min * norm_x
        elseif wᵢ > w.options.scale_max * norm_x
            wᵢ = w.options.scale_max * norm_x
        end
        w[i] = max(wᵢ, w.options.scale_abs_min)
    end
    w
end
update!(n::AbstractNorm, ::AbstractVector) = n


####################
## Specific norms ##
####################

"""
    EuclideanNorm <: AbstractNorm

The usual euclidean norm resp. 2-norm.
"""
struct EuclideanNorm <: AbstractNorm end

function distance(x::AbstractVector, y::AbstractVector, ::EuclideanNorm)
    n = length(x)
    @boundscheck n == length(y)
    @inbounds d = abs2(x[1] - y[1])
    for i in 2:n
        @inbounds d += abs2(x[i] - y[i])
    end
    sqrt(d)
end
function distance(x::AbstractVector, y::AbstractVector, w::WeightedNorm{EuclideanNorm})
    @boundscheck length(w) == length(x) == length(y)
    @inbounds d = abs2(x[1] - y[1]) / (w[1]^2)
    for i in 2:length(x)
        @inbounds d += abs2(x[i] - y[i]) / (w[i]^2)
    end
    sqrt(d)
end

function LinearAlgebra.norm(x::AbstractVector, ::EuclideanNorm)
    n = length(x)
    @inbounds d = abs2(x[1])
    for i in 2:n
        @inbounds d += abs2(x[i])
    end
    sqrt(d)
end
function LinearAlgebra.norm(x::AbstractVector, w::WeightedNorm{EuclideanNorm})
    @boundscheck length(w) == length(x)
    @inbounds out = abs2(x[1]) / (w[1]^2)
    for i in 2:length(x)
        @inbounds out += abs2(x[i]) / (w[i]^2)
    end
    sqrt(out)
end
(N::EuclideanNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::EuclideanNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    InfNorm <: AbstractNorm

The ∞-norm.
"""
struct InfNorm <: AbstractNorm end

function distance(x::AbstractVector, y::AbstractVector, ::InfNorm)
    n = length(x)
    @boundscheck n == length(y)
    @inbounds dmax = abs2(x[1] - y[1])
    for i in 2:n
        @inbounds dᵢ = abs2(x[i] - y[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
function distance(x::AbstractVector, y::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2(x[1] - y[1]) / (w[1]^2)
    for i in 2:n
        @inbounds dᵢ = abs2(x[i] - y[i]) / (w[i]^2)
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end

function LinearAlgebra.norm(x::AbstractVector, ::InfNorm)
    n = length(x)
    @inbounds dmax = abs2(x[1])
    for i in 2:n
        @inbounds dᵢ = abs2(x[i])
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
function LinearAlgebra.norm(x::AbstractVector, w::WeightedNorm{InfNorm})
    n = length(x)
    @boundscheck n == length(w)
    @inbounds dmax = abs2(x[1]) / (w[1]^2)
    for i in 2:n
        @inbounds dᵢ = abs2(x[i]) / (w[i]^2)
        dmax = Base.FastMath.max_fast(dmax, dᵢ)
    end
    sqrt(dmax)
end
(N::InfNorm)(x::AbstractVector) = LinearAlgebra.norm(x, N)
(N::InfNorm)(x::AbstractVector, y::AbstractVector) = distance(x, y, N)

"""
    euclidean_distance(u, v)

Compute ||u-v||₂.
"""
euclidean_distance(u, v) = distance(u, v, EuclideanNorm())

"""
    euclidean_norm(u)

Compute ||u||₂.
"""
euclidean_norm(x::AbstractVector) = LinearAlgebra.norm(x, EuclideanNorm())
