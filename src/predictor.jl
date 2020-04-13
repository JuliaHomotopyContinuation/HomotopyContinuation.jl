# This file contains a Padé predictor of order (2,1) used for the tracking in tracker.jl
# For a motivation of this see:
# Mixed Precision Path Tracking for Polynomial Homotopy Continuation,
# Sascha Timme (2020), arXiv:1902.02968

struct Pade21
    taylor::BitVector
    τ::Base.RefValue{Float64}
    τ_high::Base.RefValue{Float64}
    local_err::Vector{ComplexF64}
end
Pade21(n::Int) = Pade21(falses(n), Ref(Inf), Ref(Inf), zeros(ComplexF64, n))

function update!(cache::Pade21, tx::TaylorVector)
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = τ_high = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    for i in 1:length(tx)
        x, x¹, x², x³, x⁴ = tx[i]
        c¹ = abs(x¹)
        λ = max(λ_min, c¹)
        c¹ /= λ
        c² = abs(x²) / λ^2
        absx³ = abs(x³)
        c³ = absx³ / λ^3
        tol = 1e-14 * max(c¹, c², c³)
        if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
            cache.taylor[i] = true
            cache.local_err[i] = x⁴
        else
            cache.taylor[i] = false
            τᵢ = (c² / c³) / λ
            τᵢ < τ && (τ = τᵢ)
            # τ_high
            τᵢ < τ_high && (τ_high = τᵢ)
            τᵢ′ = absx³ / abs(x⁴)
            τᵢ′ < τ_high && (τ_high = τᵢ′)

            cache.local_err[i] = x⁴ - x³ * (x³ / x²)
        end
    end
    cache.τ[] = τ

    cache
end

function predict!(x̂, predictor::Pade21, H, tx::TaylorVector, Δt)
    for i in 1:length(tx)
        x, x¹, x², x³ = tx[i]
        if predictor.taylor[i]
            x̂[i] = x + Δt * (x¹ + Δt * x²)
        else
            δᵢ = 1 - Δt * x³ / x²
            x̂[i] = x + Δt * (x¹ + Δt * x² / δᵢ)
        end
    end

    nothing
end

order(c::Pade21) = 4
local_error(predictor::Pade21) = predictor.local_err
trust_region(predictor::Pade21) = predictor.τ[]
strict_trust_region(predictor::Pade21) = predictor.τ_high[]
