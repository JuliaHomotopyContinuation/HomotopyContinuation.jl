export Pade21

abstract type AbstractPredictor end
abstract type AbstractPredictorCache end

"""
    Pade21()

This uses a Padé-approximation of type (2,1) for prediction.
"""
struct Pade21 <: AbstractPredictor end

struct Pade21Cache <: AbstractPredictorCache
    x¹::Vector{ComplexF64}
    x²::Vector{ComplexF64}
    x³::Vector{ComplexF64}
    x⁴::Vector{ComplexF64}

    u::Vector{ComplexF64}
    dx¹::NTuple{1,Vector{ComplexF64}}
    dx²::NTuple{2,Vector{ComplexF64}}
    dx³::NTuple{3,Vector{ComplexF64}}

    taylor::BitVector
    τ::Base.RefValue{Float64}
    local_err::Vector{ComplexF64}
end

function cache(::Pade21, n::Int)
    x¹ = zeros(ComplexF64, n)
    x² = zeros(ComplexF64, n)
    x³ = zero(x²)
    x⁴ = zero(x²)

    u = zero(x²)
    dx¹ = (x¹,)
    dx² = (x¹, x²)
    dx³ = (x¹, x², x³)

    taylor = falses(n)
    τ = Ref(Inf)
    local_err = zero(u)

    Pade21Cache(x¹, x², x³, x⁴, u, dx¹, dx², dx³, taylor, τ, local_err)
end

function update!(cache::Pade21Cache, H, x, t, J::Jacobian, norm)
    # unpack stuff to make the rest easier to read
    @unpack u, x¹, x², x³, x⁴, dx¹, dx², dx³ = cache

    # comput all taylor coefficients x¹, x², x³, x⁴
    diff_t!(u, H, x, t)
    u .= .-u
    LA.ldiv!(x¹, J, u)

    # Check if we have to do iterative refinment for all the others as well
    δ = iterative_refinement!(x¹, J, u, norm; fixed_precision = true)
    iterative_refinement = δ > sqrt(eps())
    if iterative_refinement
        iterative_refinement!(x², J, u)
    end

    diff_t!(u, H, x, t, dx¹)
    u .= .-u
    LA.ldiv!(x², J, u)
    if iterative_refinement
        iterative_refinement!(x², J, u)
    end

    diff_t!(u, H, x, t, dx²)
    u .= .-u
    LA.ldiv!(x³, J, u)
    if iterative_refinement
        iterative_refinement!(x³, J, u)
    end

    diff_t!(u, H, x, t, dx³)
    u .= .-u
    LA.ldiv!(x⁴, J, u)
    if iterative_refinement
        iterative_refinement!(x⁴, J, u)
    end

    # This is an adaption of the algorithm outlined in
    τ = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    for i in eachindex(x)
        c¹ = fast_abs(x¹[i])
        λ = max(λ_min, c¹)
        c¹ /= λ
        c² = fast_abs(x²[i]) / λ^2
        c³ = fast_abs(x³[i]) / λ^3
        tol = 1e-14 * max(c¹, c², c³)
        if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
            cache.taylor[i] = true
            cache.local_err[i] = x⁴[i]
        else
            cache.taylor[i] = false
            τᵢ = (c² / c³) / λ
            τᵢ < τ && (τ = τᵢ)
            cache.local_err[i] = x⁴[i] - Base.FastMath.div_fast(x³[i]^2, x²[i])
        end
    end
    cache.τ[] = τ

    cache
end


function predict!(x̂, cache::Pade21Cache, H, x, t, Δt)
    @unpack x¹, x², x³, taylor = cache

    for i in eachindex(x)
        if taylor[i]
            x̂[i] = x[i] + Δt * (x¹[i] + Δt * x²[i])
        else
            δᵢ = 1 - Δt * x³[i] / x²[i]
            x̂[i] = x[i] + Δt * (x¹[i] + Δt * Base.FastMath.div_fast(x²[i], δᵢ))
        end
    end

    nothing
end

order(c::Pade21Cache) = 4
local_error(cache::Pade21Cache) = cache.local_err
trust_region(cache::Pade21Cache) = cache.τ[]
