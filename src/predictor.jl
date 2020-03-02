struct Pade21
    taylor::BitVector
    τ::Base.RefValue{Float64}
    τ_high::Base.RefValue{Float64}
    local_err::Vector{ComplexF64}
end
Pade21(n::Int) = Pade21(falses(n), Ref(Inf), Ref(Inf), zeros(ComplexF64, n))

function update!(cache::Pade21, x, x¹, x², x³, x⁴)
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = τ_high = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    for i in eachindex(x)
        c¹ = abs(x¹[i])
        λ = max(λ_min, c¹)
        c¹ /= λ
        c² = abs(x²[i]) / λ^2
        absx³ = abs(x³[i])
        c³ = absx³ / λ^3
        tol = 1e-14 * max(c¹, c², c³)
        if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
            cache.taylor[i] = true
            cache.local_err[i] = x⁴[i]
        else
            cache.taylor[i] = false
            τᵢ = (c² / c³) / λ
            τᵢ < τ && (τ = τᵢ)
            # τ_high
            τᵢ < τ_high && (τ_high = τᵢ)
            τᵢ′ = absx³ / abs(x⁴[i])
            τᵢ′ < τ_high && (τ_high = τᵢ′)

            cache.local_err[i] = x⁴[i] - x³[i] * (x³[i] / x²[i])
        end
    end
    cache.τ[] = τ

    cache
end

function predict!(x̂, predictor::Pade21, H, x, x¹, x², x³, Δt)
    for i in eachindex(x)
        if predictor.taylor[i]
            x̂[i] = x[i] + Δt * (x¹[i] + Δt * x²[i])
        else
            δᵢ = 1 - Δt * x³[i] / x²[i]
            x̂[i] = x[i] + Δt * (x¹[i] + Δt * x²[i] / δᵢ)
        end
    end

    nothing
end

order(c::Pade21) = 4
local_error(predictor::Pade21) = predictor.local_err
trust_region(predictor::Pade21) = predictor.τ[]
strict_trust_region(predictor::Pade21) = predictor.τ_high[]
