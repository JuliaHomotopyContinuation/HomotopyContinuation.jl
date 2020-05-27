export Pade21

"""
    Pade21()

This uses a Padé-approximation of type (2,1) for prediction.
"""
struct Pade21 <: AbstractPredictor end
struct Pade21Cache{T,AV<:AbstractVector{T}} <: AbstractPredictorCache
    x²::Vector{T}
    x³::Vector{T}
    x_h::AV

    u::Vector{T}
    u₁::Vector{T}
    u₂::Vector{T}
end

function cache(::Pade21, H, x, ẋ, t)
    x², x³, x_h = copy(ẋ), copy(ẋ), copy(x)
    u = evaluate(H, x, t)
    u₁, u₂ = copy(u), copy(u)
    Pade21Cache(x², x³, x_h, u, u₁, u₂)
end

@inline function g!(u, H, x, ẋ, t, h, x_h)
    @inbounds for i in eachindex(x)
        x_h[i] = x[i] + h * ẋ[i]
    end
    evaluate!(u, H, x_h, t + h)
end

@inline function g!(u, H, x, ẋ, ẍ, t, h, x_h)
    h² = 0.5 * h * h
    @inbounds for i in eachindex(x)
        x_h[i] = x[i] + h * ẋ[i] + h² * ẍ[i]
    end
    evaluate!(u, H, x_h, t + h)
end

function update!(cache::Pade21Cache, H, x, ẋ, t, Jac::JacobianMonitor, eval_err::Float64)
    # unpack stuff to make the rest easier to read
    @unpack u, u₁, u₂, x_h, x², x³ = cache
    ψ = 10eval_err
    h₂ = nthroot(ψ, 4)
    g!(u₁, H, x, ẋ, t, h₂, x_h)
    g!(u₂, H, x, ẋ, t, -h₂, x_h)

    h² = h₂^2
    @inbounds for i in eachindex(u₁)
        u[i] = -(u₁[i] + u₂[i]) / h²
    end
    LA.ldiv!(x², Jac, u)

    h₃ = nthroot(ψ, 5)
    g!(u₁, H, x, ẋ, x², t, h₃, x_h)
    g!(u₂, H, x, ẋ, x², t, -h₃, x_h)

    h³ = h₃^3
    @inbounds for i in eachindex(u₁)
        u[i] = -3 * (u₁[i] - u₂[i]) / h³
    end
    LA.ldiv!(x³, Jac, u)

    cache
end


function predict!(x̂, cache::Pade21Cache, H::HomotopyWithCache, x, t, Δt, ẋ, J)
    @unpack x², x³ = cache
    @inbounds for i in eachindex(x)
        δ = 1 - Δt * x³[i] / (3 * x²[i])
        if isnan(δ) || iszero(δ)
            x̂[i] = x[i] + Δt * ẋ[i] + 0.5 * Δt^2 * x²[i]
        else
            x̂[i] = x[i] + Δt * ẋ[i] + (0.5 * Δt^2 * x²[i]) / δ
        end
    end
    nothing
end

order(::Pade21Cache) = 4
@inline highest_derivative(cache::Pade21Cache) = (cache.x³, 3)
second_derivative(cache::Pade21Cache) = cache.x²
third_derivative(cache::Pade21Cache) = cache.x³
