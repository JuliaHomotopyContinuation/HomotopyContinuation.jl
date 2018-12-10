import LinearAlgebra
using ..Utilities

export Newton

"""
    Newton()

A classical simple Newton operator for square linear systems using the LU factorization
to solve the linear systems.
"""
struct Newton <: AbstractCorrector end

struct NewtonCache{T} <: AbstractCorrectorCache
    Jᵢ::Matrix{T}
    rᵢ::Vector{T}
end

function cache(::Newton, H::HomotopyWithCache, x, t)
    Jᵢ = Homotopies.jacobian(H, x, t)
    rᵢ = Homotopies.evaluate(H, x, t)
    NewtonCache(Jᵢ, rᵢ)
end


function correct!(out, alg::Newton, cache::NewtonCache, H::HomotopyWithCache, x₀, t, tol, maxit)
    Jᵢ, rᵢ= cache.Jᵢ, cache.rᵢ

    copyto!(out, x₀)
    xᵢ₊₁ = xᵢ = out # just alias to make logic easier
    rᵢ₊₁ = rᵢ

    # initialize variables
    T = real(eltype(xᵢ))
    Θ₀ = Θᵢ₋₁ = norm_Δxᵢ₋₁ = norm_Δxᵢ = norm_Δx₀ = zero(T)
    accuracy = T(Inf)
    ω₀ = ω = 0.0
    for i ∈ 0:(maxit-1)
        evaluate_and_jacobian!(rᵢ, Jᵢ, H, xᵢ, t)
        Δxᵢ = solve!(Jᵢ, rᵢ)
        norm_Δxᵢ₋₁ = norm_Δxᵢ
        accuracy = norm_Δxᵢ = euclidean_norm(Δxᵢ)
        for k in eachindex(xᵢ)
            xᵢ₊₁[k] = xᵢ[k] - Δxᵢ[k]
        end

        if i == 0
            norm_Δx₀ = norm_Δxᵢ₋₁ = norm_Δxᵢ
            if norm_Δx₀ ≤ tol
                return Result(converged, norm_Δx₀, i + 1, 0.0, 0.0, norm_Δx₀)
            end
            continue
        end

        ωᵢ₋₁ = 2norm_Δxᵢ / norm_Δxᵢ₋₁^2
        if i == 1
            ω = ω₀ = ωᵢ₋₁
        else
            ω = @fastmath max(ω, ωᵢ₋₁)
        end

        if accuracy ≤ tol
            return Result(converged, accuracy, i + 1, ω₀, ω, norm_Δx₀)
        end
    end

    return Result(maximal_iterations, accuracy, maxit, ω₀, ω, norm_Δx₀)
end
