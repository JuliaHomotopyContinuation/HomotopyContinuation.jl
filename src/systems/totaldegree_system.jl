export TotalDegreeSystem

"""
    TotalDegreeSystem(polynomials, vars) <: AbstractSystem

Create a tottal degree system
"""
struct TotalDegreeSystem{T} <: AbstractSystem
    degrees::Vector{Int}
    degree_idxs::Vector{Int}
    hom_idx::Int # 0 if affine
    scaling_factors::Vector{T}

    function TotalDegreeSystem(
        degrees::Vector{Int},
        degree_idxs::Vector{Int},
        hom_idx::Int = 0,
        scaling_factors::Vector{T} = ones(length(degrees));
        affine = false,
    ) where {T}
        @assert length(degrees) == length(degree_idxs) == length(scaling_factors)
        if affine
            hom_idx = 0
        else
            @assert 1 ≤ hom_idx ≤ length(degrees) + 1
        end
        new{T}(degrees, degree_idxs, hom_idx, scaling_factors)
    end
end
function TotalDegreeSystem(F::Vector{<:MP.AbstractPolynomialLike}, args...; kwargs...)
    TotalDegreeSystem(MP.maxdegree.(F), args...; kwargs...)
end
function TotalDegreeSystem(degrees::Vector{Int}, args...; kwargs...)
    TotalDegreeSystem(
        degrees,
        collect(1:length(degrees)),
        length(degrees) + 1,
        args...;
        kwargs...,
    )
end

function Base.size(F::TotalDegreeSystem)
    n = length(F.degrees)
    F.hom_idx == 0 ? (n, n) : (n, n + 1)
end

cache(::TotalDegreeSystem, x) = SystemNullCache()

function evaluate!(u, F::TotalDegreeSystem, x::ProjectiveVectors.PVector, ::SystemNullCache)
    λ = F.scaling_factors
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        u[i] = λ[i] * (x[i]^d - x[end]^d)
    end
    u
end

function evaluate!(u, F::TotalDegreeSystem, x::AbstractVector, ::SystemNullCache)
    λ = F.scaling_factors
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        dix = F.degree_idxs[i]
        u[i] = λ[i] * x[dix]^d - λ[i]
    end
    u
end
function evaluate(F::TotalDegreeSystem, x, cache::SystemNullCache)
    u = similar(x, size(F, 1))
    evaluate!(u, F, x, cache)
    u
end

function jacobian!(U, F::TotalDegreeSystem, x::PVector{<:Number,1}, ::SystemNullCache)
    U .= zero(eltype(x))
    λ = F.scaling_factors
    N = length(x)
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            U[i, didx] = λ[i]
            U[i, N] = -λ[i]
        elseif d > 1
            U[i, didx] = λ[i] * d * x[didx]^(d - 1)
            U[i, N] = -λ[i] * d * x[N]^(d - 1)
        end
    end
    U
end
function jacobian!(U, F::TotalDegreeSystem, x::AbstractVector, ::SystemNullCache)
    U .= zero(eltype(x))
    λ = F.scaling_factors
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            U[i, didx] = λ[i]
        elseif d > 1
            U[i, didx] = λ[i] * d * x[didx]^(d - 1)
        end
    end
    U
end

function jacobian(F::TotalDegreeSystem, x, cache::SystemNullCache)
    U = similar(x, size(F))
    jacobian!(U, F, x, cache)
    U
end

function evaluate_and_jacobian!(
    u,
    U,
    F::TotalDegreeSystem,
    x::PVector{<:Number,1},
    ::SystemNullCache,
)
    U .= zero(eltype(x))
    N = length(x)
    λ = F.scaling_factors
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            u[i] = λ[i] * (x[didx] - x[N])
            U[i, didx] = λ[i]
            U[i, N] = -λ[i]
        elseif d > 1
            xd = x[didx]^(d - 1)
            xh = x[N]^(d - 1)
            u[i] = λ[i] * (xd * x[didx] - xh * x[N])
            U[i, didx] = λ[i] * d * xd
            U[i, N] = -λ[i] * d * xh
        end
    end
    nothing
end
function evaluate_and_jacobian!(
    u,
    U,
    F::TotalDegreeSystem,
    x::AbstractVector,
    ::SystemNullCache,
)
    U .= zero(eltype(x))
    λ = F.scaling_factors
    @inbounds for i = 1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            u[i] = λ[i] * x[didx] - λ[i]
            U[i, didx] = λ[i]
        elseif d > 1
            xd = x[didx]^(d - 1)
            u[i] = λ[i] * xd * x[didx] - λ[i]
            U[i, didx] = λ[i] * d * xd
        end
    end
    nothing
end

function evaluate_and_jacobian(F::TotalDegreeSystem, x, cache::SystemNullCache)
    u = similar(x, size(F, 1))
    U = similar(x, size(F))
    evaluate_and_jacobian!(u, U, F, x, cache)
    u, U
end
