import LinearAlgebra
import Random
using ..Utilities

export Newton

"""
    Newton(;simplified_last_step=true)

An ordinary Newton's method. If `simplified_last_step` is `true`, then for the last iteration
the previously Jacobian will be used. This uses an LU-factorization for square systems
and a QR-factorization for overdetermined.
"""
struct Newton <: AbstractCorrector
    simplified_last_step::Bool
end

Newton(;simplified_last_step=true) = Newton(simplified_last_step)

struct NewtonCache{T, Fac<:LinearAlgebra.Factorization} <: AbstractCorrectorCache
    Jac::Jacobian{T, Fac}
    rᵢ::Vector{T}
    Δxᵢ::Vector{T}
end

function cache(::Newton, H::HomotopyWithCache, x, t)
    Jac = Jacobian(Homotopies.jacobian(H, x, t))
    rᵢ = Homotopies.evaluate(H, x, t)
    Δxᵢ = copy(rᵢ)

    NewtonCache(Jac, rᵢ, Δxᵢ)
end

function correct!(out, alg::Newton, cache::NewtonCache, H::HomotopyWithCache, x₀, t; tol=1e-6, maxiters::Integer=3)
    Jac, rᵢ, Δxᵢ = cache.Jac, cache.rᵢ, cache.Δxᵢ
    Jᵢ = Jac.J
    copyto!(out, x₀)
    xᵢ₊₁ = xᵢ = out # just alias to make logic easier
    rᵢ₊₁ = rᵢ

    # initialize variables
    T = real(eltype(xᵢ))
    Θ₀ = Θᵢ₋₁ = norm_Δxᵢ₋₁ = norm_Δxᵢ = norm_Δx₀ = zero(T)
    accuracy = T(Inf)
    ω₀ = ω = 0.0
    for i ∈ 0:(maxit)
        if i == maxit && alg.simplified_last_step
            evaluate!(rᵢ, H, xᵢ, t)
        else
            evaluate_and_jacobian!(rᵢ, Jᵢ, H, xᵢ, t)
            Utilities.updated_jacobian!(Jac)
        end
        solve!(Δxᵢ, Jac, rᵢ)
        norm_Δxᵢ₋₁ = norm_Δxᵢ
        norm_Δxᵢ = euclidean_norm(Δxᵢ)
        @inbounds for k in eachindex(xᵢ)
            xᵢ₊₁[k] = xᵢ[k] - Δxᵢ[k]
        end

        if i == 0
            norm_Δx₀ = norm_Δxᵢ₋₁ = norm_Δxᵢ
            if norm_Δx₀ ≤ tol
                return Result(converged, norm_Δx₀, i + 1, 0.0, 0.0, norm_Δx₀)
            end
            continue
        end

        Θᵢ₋₁ = norm_Δxᵢ / norm_Δxᵢ₋₁
        ωᵢ₋₁ = 2Θᵢ₋₁ / norm_Δxᵢ₋₁
        if i == 1
            ω = ω₀ = ωᵢ₋₁
        else
            ω = @fastmath max(ω, ωᵢ₋₁)
        end
        if Θᵢ₋₁ > 0.5
            return Result(terminated, norm_Δx₀, i + 1, ω₀, ω, norm_Δx₀)
        end

        accuracy = norm_Δxᵢ / (1 - 2Θᵢ₋₁^2)
        if accuracy ≤ tol
            return Result(converged, accuracy, i + 1, ω₀, ω, norm_Δx₀)
        end
    end

    return Result(maximal_iterations, accuracy, maxiters, ω₀, ω, norm_Δx₀)
end
