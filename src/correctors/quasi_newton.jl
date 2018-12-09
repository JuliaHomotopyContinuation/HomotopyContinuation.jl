export QuasiNewton

import ..Utilities
import LinearAlgebra

"""
    QuasiNewton()

A classical simple Newton operator for square linear systems using the LU factorization
to solve the linear systems.
"""
struct QuasiNewton <: AbstractCorrector
    k_max::Int
end
QuasiNewton(;k_max=100) = QuasiNewton(k_max)

struct QuasiNewtonCache{T, R<:Real} <: AbstractCorrectorCache
    Jᵢ::Matrix{T}
    rᵢ::Vector{T}
    factorization::LinearAlgebra.LU{T, Matrix{T}}
    σ::Vector{R}
    δx::Matrix{T}
end

function cache(alg::QuasiNewton, H::HomotopyWithCache, x, t)
    Jᵢ = Homotopies.jacobian(H, x, t)
    rᵢ = Homotopies.evaluate(H, x, t)
    factorization = LinearAlgebra.generic_lufact!(Jᵢ)
    σ = zeros(typeof(abs2(rᵢ[1])), alg.k_max)
    δx = similar(Jᵢ, length(rᵢ), alg.k_max)
    QuasiNewtonCache(Jᵢ, rᵢ, factorization, σ, δx)
end

function correct!(out, alg::QuasiNewton, cache::QuasiNewtonCache, H::HomotopyWithCache, x₀, t, tol, maxit, patch_state)
    Jᵢ, rᵢ, σ, δx = cache.Jᵢ, cache.rᵢ, cache.σ, cache.δx
    LU = cache.factorization
    @inbounds begin
        for i in eachindex(x₀)
            out[i] = x₀[i]
        end

        # initialize variables
        xᵢ₊₁ = xᵢ = out # just assign to make logic easier
        Θ₀ = Θᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = zero(real(eltype(xᵢ)))
        accuracy = convert(typeof(norm_Δxᵢ₋₁), Inf)

        evaluate_and_jacobian!(rᵢ, Jᵢ, H, xᵢ, t)
        for i in eachindex(rᵢ)
            rᵢ[i] = -rᵢ[i]
        end
        # we only need to factorize once
        Utilities.fast_factorization!(LU)

        δxᵢ₊₁ = δxᵢ = δx₀ = Utilities.fast_ldiv!(LU, rᵢ)
        colset!(δx, δx₀, 1)
        σᵢ₊₁ = σᵢ = σ₀ = σ[1] = fast_2_norm2(δx₀)
        for i in eachindex(xᵢ₊₁)
            xᵢ₊₁[i] = xᵢ[i] + δxᵢ[i]
        end
        norm_Δx₀ = sqrt(σ₀)
        if norm_Δx₀ < tol
            return Result(converged, norm_Δx₀, 1, Θ₀, Θᵢ, norm_Δx₀)
        end

        for k ∈ 0:maxit
            rᵢ₊₁ = evaluate!(rᵢ, H, xᵢ₊₁, t)
            for i in eachindex(rᵢ₊₁)
                rᵢ₊₁[i] = -rᵢ₊₁[i]
            end

            v = Utilities.fast_ldiv!(LU, rᵢ₊₁)
            for i ∈ 1:k
                ᾱ = coldot(δx, i, v) / σ[i]
                colmuladd!(v, ᾱ, δx, i + 1)
            end
            αᵢ₊₁ = coldot(δx, k+1, v) / σ[k+1]
            Θᵢ = sqrt(fast_2_norm2(v) / σ[k+1])
            if k == 0
                Θ₀ = Θᵢ
            end
            if Θᵢ ≥ 0.5
                return Result(terminated, accuracy, k + 1, Θ₀, Θᵢ, norm_Δx₀)
            end
            λ = (@fastmath 1.0 / (1.0 - αᵢ₊₁))
            for i in eachindex(δxᵢ₊₁)
                δxᵢ₊₁[i] = λ * v[i]
            end
            colset!(δx, δxᵢ₊₁, k+2)
            σᵢ₊₁ = σ[k+2] = fast_2_norm2(δxᵢ₊₁)
            for i in eachindex(xᵢ₊₁)
                xᵢ₊₁[i] += δxᵢ₊₁[i]
            end
            accuracy = sqrt(σᵢ₊₁)
            if accuracy ≤ tol
                return Result(converged, accuracy, k + 1, Θ₀, Θᵢ, norm_Δx₀)
            end
        end
        return Result(maximal_iterations, accuracy, maxit, Θ₀, Θᵢ, norm_Δx₀)
    end
end

"""
    colset!(A, v, j)

Computes `A[:,j] .= v`.
"""
Base.@propagate_inbounds function colset!(A::AbstractMatrix, v::AbstractVector, j::Integer)
    @boundscheck length(v) == size(A, 1) && 1 ≤ j ≤ size(A, 2)

    @inbounds for i ∈ eachindex(v)
        A[i, j] = v[i]
    end
end


"""
    coldot(v, A, j)

Computes `dot(v, A[:, j])`.
"""
Base.@propagate_inbounds function coldot(v::AbstractVector, A::AbstractMatrix, j::Integer)
    @boundscheck length(v) == size(A, 1) && 1 ≤ j ≤ size(A, 2)

    out = zero(eltype(v))
    @inbounds for i ∈ eachindex(v)
        out += conj(v[i]) * A[i, j]
    end
    out
end

"""
    coldot(A, j, v)

Computes `dot(A[:, j], v)`.
"""
Base.@propagate_inbounds function coldot(A::AbstractMatrix, j::Integer, v::AbstractVector)
    @boundscheck length(v) == size(A, 1) && 1 ≤ j ≤ size(A, 2)

    out = zero(eltype(v))
    @inbounds for i ∈ eachindex(v)
        out += conj(A[i, j]) * v[i]
    end
    out
end

"""
    colmuladd!(v, α, A, j)

Computes `v .= v .+ α .* A[:, j]`.
"""
Base.@propagate_inbounds function colmuladd!(v::AbstractVector, α, A::AbstractMatrix, j::Integer)
    @boundscheck length(v) == size(A, 1) && 1 ≤ j ≤ size(A, 2)

    @inbounds for i ∈ eachindex(v)
        v[i] += α * A[i, j]
    end
    v
end
