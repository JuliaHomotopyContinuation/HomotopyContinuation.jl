import LinearAlgebra

export Taylor

"""
    Taylor()
This uses the explicit Taylor method for prediction, also known as the
tangent predictor.
"""
struct Taylor <: AbstractStatelessPredictor end
struct TaylorCache{T} <: AbstractStatelessPredictorCache
    x²::Vector{T}
    x³::Vector{T}
    x_h::Vector{T}
    u₁::Vector{T}
    u₂::Vector{T}
    u₃::Vector{T}
    u₄::Vector{T}

    h₂::Float64
    h₃::Float64
end

function cache(::Taylor, H, x, t)
    x², x³, x_h = copy(x.data), copy(x.data), copy(x.data)
    u₁, u₂, u₃, u₄ = copy(x.data), copy(x.data), copy(x.data), copy(x.data)
    h₂ = nthroot(eps(), 2+4) # x² is the second order derivative and has a 4th order approximation
    h₃ = nthroot(eps(), 3+2) # x³ is the third order derivative and has a 2th order approximation
    TaylorCache(x², x³, x_h, u₁, u₂, u₃, u₄, h₂, h₃)
end

function g₁!(u, H, x, ẋ, t, h, x_h)
    @inbounds for i in eachindex(x)
        x_h[i] = x[i] + h * ẋ[i]
    end
    Homotopies.evaluate!(u, H, x_h, t + h)
end

function g₂!(u, H, x, ẋ, x², t, h, x_h)
    h² = 0.5 * h * h
    @inbounds for i in eachindex(x)
        x_h[i] = x[i] + h * ẋ[i] + h² * x²[i]
    end
    Homotopies.evaluate!(u, H, x_h, t + h)
end

function update!(cache::TaylorCache, H, x, ẋ, t, fac)
    # unpack stuff to make the rest easier to read
    u₁, u₂, u₃, u₄ = cache.u₁, cache.u₂, cache.u₃, cache.u₄
    x_h, h₂, h₃ = cache.x_h, cache.h₂, cache.h₃
    x², x³ = cache.x², cache.x³

    # compute the second derivative at (x,t)
    g₁!(u₁, H, x, ẋ, t, h₂, x_h)
    g₁!(u₂, H, x, ẋ, t, -h₂, x_h)
    g₁!(u₃, H, x, ẋ, t, 2h₂, x_h)
    g₁!(u₄, H, x, ẋ, t, -2h₂, x_h)

    h² = h₂^2
    @inbounds for i in eachindex(u₁)
        x²[i] = (-4/3 * (u₁[i] + u₂[i]) + (u₃[i] + u₄[i]) / 12) / h²
    end
    solve!(fac, x²)


    g₂!(u₁, H, x, ẋ, x², t, h₃, x_h)
    g₂!(u₂, H, x, ẋ, x², t, -h₃, x_h)

    h³ = h₃ * h₃ * h₃
    @inbounds for i in eachindex(u₁)
        x³[i] = -3(u₁[i] - u₂[i]) / h³
    end
    solve!(fac, x³)

    cache
end

function predict!(xnext, ::Taylor, cache::TaylorCache, H::HomotopyWithCache, x, t, Δt, ẋ)
    x², x³ = cache.x², cache.x³
    δ = cache.x_h
    @inbounds for i in eachindex(x)
        δ[i] = -x³[i] / (3 * x²[i])
    end

    @inbounds for i in eachindex(x)
        xnext[i] = x[i] + Δt * ẋ[i] + (0.5 * Δt^2 * x²[i]) / (1 + δ[i]*Δt)
    end
    nothing
end

order(::Taylor) = 4
