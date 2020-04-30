# This file contains a Padé predictor of order (2,1) used for the tracking in tracker.jl
# For a motivation of this see:
# Mixed Precision Path Tracking for Polynomial Homotopy Continuation,
# Sascha Timme (2020), arXiv:1902.02968

mutable struct Predictor
    taylor::BitVector
    τ::Float64
    local_err::Vector{ComplexF64}
    ord::Int
end
Predictor(n::Int) = Predictor(falses(n), Inf, zeros(ComplexF64, n), 4)

function update!(pred::Predictor, tx::TaylorVector, trust_tx::Vector{Bool})
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    if trust_tx[1] && trust_tx[2] && trust_tx[3]
        for i = 1:length(tx)
            x, x¹, x², x³, x⁴ = tx[i]
            c¹ = abs(x¹)
            λ = max(λ_min, c¹)
            c¹ /= λ
            c² = abs(x²) / λ^2
            absx³ = abs(x³)
            c³ = absx³ / λ^3
            tol = 1e-14 * max(c¹, c², c³)
            if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
                pred.taylor[i] = true
                pred.local_err[i] = x⁴
            else
                pred.taylor[i] = false
                τᵢ = (c² / c³) / λ
                τᵢ < τ && (τ = τᵢ)
                pred.local_err[i] = x⁴ - x³ * (x³ / x²)
            end
        end
        pred.ord = 4
        pred.τ = τ
    else
        _, x¹, x², x³, = vectors(tx)
        pred.taylor .= true
        if trust_tx[2]
            pred.ord = 3
            pred.local_err .= x³
        else
            pred.ord = 2
            pred.local_err .= x²
        end
    end

    pred.τ = τ

    pred
end

function predict!(x̂, predictor::Predictor, H, tx::TaylorVector, Δt)
    for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        if predictor.taylor[i]
            if predictor.ord ≥ 3
                x̂[i] = x + Δt * (x¹ + Δt * x²)
            else
                x̂[i] = x + Δt * x¹
            end
        else
            δᵢ = 1 - Δt * x³ / x²
            x̂[i] = x + Δt * (x¹ + Δt * x² / δᵢ)
        end
    end

    nothing
end

order(P::Predictor) = P.ord
local_error(predictor::Predictor) = predictor.local_err
trust_region(predictor::Predictor) = predictor.τ
