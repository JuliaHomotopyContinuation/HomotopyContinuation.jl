export Pade21

"""
    Pade21()

This uses a Padé-approximation of type (2,1) for prediction.
"""
struct Pade21 <: AbstractPredictor end
struct Pade21Cache{T, AV<:AbstractVector{T}} <: AbstractPredictorCache
    x²::Vector{T}
    x³::Vector{T}
    x_h::AV

    u::Vector{T}
    u₁::Vector{T}
    u₂::Vector{T}

    h₂::Float64
    h₃::Float64
end

function cache(::Pade21, H, x, ẋ, t)
    x², x³, x_h = copy(ẋ), copy(ẋ), copy(x)
    u = evaluate(H, x, t)
    u₁, u₂ = copy(u), copy(u)
    h₂ = nthroot(eps(), 2+2) # x² is the second order derivative and has a 2nd order approximation
    h₃ = nthroot(eps(), 3+2) # x³ is the third order derivative and has a 2nd order approximation
    Pade21Cache(x², x³, x_h, u, u₁, u₂, h₂, h₃)
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

function update!(cache::Pade21Cache, H, x, ẋ, t, Jac::Jacobian)
    # unpack stuff to make the rest easier to read
    u, u₁, u₂ = cache.u, cache.u₁, cache.u₂
    x_h, h₂, h₃ = cache.x_h, cache.h₂, cache.h₃
    x², x³ = cache.x², cache.x³

    h₂ *= max(abs(t), 1.0)
    # compute the second derivative at (x,t)
    g!(u₁, H, x, ẋ, t, h₂, x_h)
    g!(u₂, H, x, ẋ, t, -h₂, x_h)

    h² = h₂^2
    @inbounds for i in eachindex(u₁)
        u[i] = -(u₁[i] + u₂[i]) / h²
    end
    solve!(x², Jac, u)

    h₃ *= max(abs(t), 1.0)
    g!(u₁, H, x, ẋ, x², t, h₃, x_h)
    g!(u₂, H, x, ẋ, x², t, -h₃, x_h)

    h³ = h₃^3
    @inbounds for i in eachindex(u₁)
        u[i] = -3(u₁[i] - u₂[i]) / h³
    end
    solve!(x³, Jac, u)

    cache
end

function predict!(xnext, ::Pade21, cache::Pade21Cache, H::HomotopyWithCache, x, t, Δt, ẋ, Jac::Jacobian)
    x², x³ = cache.x², cache.x³
    @inbounds for i in eachindex(x)
        δ = 1 - Δt * x³[i] / (3 * x²[i])
        if isnan(δ) || iszero(δ)
            xnext[i] = x[i] + Δt * ẋ[i] + 0.5 * Δt^2 * x²[i]
        else
            xnext[i] = x[i] + Δt * ẋ[i] + (0.5 * Δt^2 * x²[i]) / δ
        end
    end
    nothing
end

order(::Pade21) = 4
