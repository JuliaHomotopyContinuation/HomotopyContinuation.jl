import LinearAlgebra

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
    Δx̄ᵢ₊₁::Vector{T}
    factorization::LinearAlgebra.LU{T, Matrix{T}}
end

function cache(::Newton, H::HomotopyWithCache, x, t)
    Jᵢ = Homotopies.jacobian(H, x, t)
    rᵢ = Homotopies.evaluate(H, x, t)
    Δx̄ᵢ₊₁ = copy(rᵢ)
    factorization = LinearAlgebra.LU(copy(Jᵢ), ones(Int, size(Jᵢ, 1)), 0)
    NewtonCache(Jᵢ, rᵢ, Δx̄ᵢ₊₁, factorization)
end


function correct!(out, alg::Newton, cache::NewtonCache, H::HomotopyWithCache, x₀, t, tol, maxit, patch_state)
    Jᵢ, rᵢ, factorization, Δx̄ᵢ₊₁ = cache.Jᵢ, cache.rᵢ, cache.factorization, cache.Δx̄ᵢ₊₁

    copyto!(out.data, x₀.data)
    xᵢ₊₁ = xᵢ = out # just assign to make logic easier
    rᵢ₊₁ = rᵢ

    # initialize variables
    T = real(eltype(xᵢ))
    Θ₀ = Θᵢ = norm_Δxᵢ = norm_Δx₀ = zero(T)
    accuracy = T(Inf)

    @inbounds for i ∈ 0:(maxit-1)
        if i > 0
            AffinePatches.update!(patch_state, xᵢ, true)
        end
        if i == 0
            evaluate_and_jacobian!(rᵢ, Jᵢ, H, xᵢ, t)
        else
            jacobian!(Jᵢ, H, xᵢ, t)
        end
        fast_factorization!(factorization, Jᵢ)
        Δxᵢ = fast_ldiv!(factorization, rᵢ)

        norm_Δxᵢ = fast_2_norm(Δxᵢ)
        # Update
        for k in eachindex(xᵢ)
            xᵢ₊₁[k] = xᵢ[k] - Δxᵢ[k]
        end
        AffinePatches.normalize!(xᵢ₊₁, patch_state)

        evaluate!(rᵢ₊₁, H, xᵢ₊₁, t)
        copyto!(Δx̄ᵢ₊₁, rᵢ₊₁)
        fast_ldiv!(factorization, Δx̄ᵢ₊₁)
        accuracy = norm_Δx̄ᵢ₊₁ = fast_2_norm(Δx̄ᵢ₊₁)
        Θᵢ = norm_Δx̄ᵢ₊₁ / norm_Δxᵢ
        if i == 0
            Θ₀ = Θᵢ
            norm_Δx₀ = norm_Δxᵢ
        end

        if accuracy ≤ tol || (i == 0 && norm_Δxᵢ ≤ tol)
            return Result(converged, accuracy, i + 1, Θ₀, Θᵢ, norm_Δx₀)
        end

        if Θᵢ ≥ 1.0
            return Result(terminated, accuracy, i + 1, Θ₀, Θᵢ, norm_Δx₀)
        end
    end

    return Result(maximal_iterations, accuracy, maxit, Θ₀, Θᵢ, norm_Δx₀)
end
