export TotalDegreeSystem

"""
    TotalDegreeSystem(polynomials, vars) <: AbstractSystem

Create a tottal degree system
"""
struct TotalDegreeSystem{T} <: AbstractSystem
    degrees::Vector{Int}
    degree_idxs::Vector{Int}
    hom_idx::Int
    scaling_factors::Vector{T}

    function TotalDegreeSystem(degrees::Vector{Int}, degree_idxs::Vector{Int}, hom_idx::Int,
            scaling_factors::Vector{T}=ones(length(degrees))) where {T}
        @assert length(degrees) == length(degree_idxs) == length(scaling_factors)
        @assert 1 ≤ hom_idx ≤ length(degrees) + 1
        new{T}(degrees, degree_idxs, hom_idx, scaling_factors)
    end
end
function TotalDegreeSystem(F::Vector{<:MP.AbstractPolynomialLike}, args...)
    TotalDegreeSystem(MP.maxdegree.(F), args...)
end
function TotalDegreeSystem(degrees::Vector{Int}, args...)
    TotalDegreeSystem(degrees, collect(1:length(degrees)), length(degrees) + 1, args...)
end

Base.size(F::TotalDegreeSystem) = (length(F.degrees), length(F.degrees) + 1)

cache(::TotalDegreeSystem, x) = SystemNullCache()

function evaluate!(u, F::TotalDegreeSystem, x, ::SystemNullCache)
    λ = F.scaling_factors
    @inbounds for i=1:length(F.degrees)
        d = F.degrees[i]
        dix = F.degree_idxs[i]
        u[i] = λ[i] * (x[dix]^d - x[F.hom_idx]^d)
    end
    u
end
function evaluate(F::TotalDegreeSystem, x, cache::SystemNullCache)
    u = similar(x, size(F, 1))
    evaluate!(u, F, x, cache)
    u
end

function jacobian!(U, F::TotalDegreeSystem, x, ::SystemNullCache)
    U .= zero(eltype(x))
    hidx = F.hom_idx
    λ = F.scaling_factors
    @inbounds for i=1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            U[i, didx] = λ[i]
            U[i, hidx] = -λ[i]
        elseif d > 1
            U[i, didx] = λ[i] * d * x[didx]^(d-1)
            U[i, hidx] = -λ[i] * d * x[hidx]^(d-1)
        end
    end
    U
end
function jacobian(F::TotalDegreeSystem, x, cache::SystemNullCache)
    U = similar(x, size(F))
    jacobian!(U, F, x, cache)
    U
end

function evaluate_and_jacobian!(u, U, F::TotalDegreeSystem, x, ::SystemNullCache)
    U .= zero(eltype(x))
    hidx = F.hom_idx
    λ = F.scaling_factors
    @inbounds for i=1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            u[i] = λ[i] * (x[didx] - x[hidx])
            U[i, didx] = λ[i]
            U[i, hidx] = -λ[i]
        elseif d > 1
            xd = x[didx]^(d-1)
            xh = x[hidx]^(d-1)
            u[i] = λ[i] * (xd * x[didx] - xh * x[hidx])
            U[i, didx] = λ[i] * d * xd
            U[i, hidx] = -λ[i] * d * xh
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


weylnorm2(F::TotalDegreeSystem) = 2 * length(F.degrees)
