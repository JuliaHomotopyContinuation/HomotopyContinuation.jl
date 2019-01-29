export FixedPointHomotopy, FixedPointHomotopyCache, gamma

"""
    FixedPointHomotopy(F, x₀; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = (1-t)F(x) + γt(x-x₀)``.
"""
struct FixedPointHomotopy{S<:AbstractSystem,T<:AbstractVector} <: AbstractHomotopy
    F::S
    x₀::T
    gamma::Complex{Float64}

end
function FixedPointHomotopy(F::AbstractSystem, x₀::AbstractVector; gamma=randomish_gamma())
    FixedPointHomotopy(F, x₀, gamma)
end

"""
    FixedPointHomotopyCache

An simple cache for `StartTargetHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct FixedPointHomotopyCache{SC} <: AbstractHomotopyCache
    F::SC
end

function cache(H::FixedPointHomotopy, x, t)
    F_cache = cache(H.F, x)
    FixedPointHomotopyCache(F_cache)
end

Base.size(H::FixedPointHomotopy) = size(H.F)

"""
    gamma(H::FixedPointHomotopy)

Obtain the gamma used in the FixedPointHomotopy.
"""
gamma(H::FixedPointHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the FixedPointHomotopy.
"""
γ(H::FixedPointHomotopy) = gamma(H)

function evaluate!(u, H::FixedPointHomotopy, x, t, c::FixedPointHomotopyCache)
    evaluate!(u, H.F, x, c.F)

    u .= (γ(H) * t) .* (x .- H.x₀) .+ (1 - t) .* u

    u
end
(H::FixedPointHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function dt!(u, H::FixedPointHomotopy, x, t, c::FixedPointHomotopyCache)
    evaluate!(u, H.F, x, c.F)
    u .= γ(H) .* (x .- H.x₀) .- u
    u
end

function jacobian!(U, H::FixedPointHomotopy, x, t, c::FixedPointHomotopyCache)
    jacobian!(U, H.F, x, c.F)

    U .= (1 - t) .* U
    γt = γ(H) * t
    for i=1:length(x)
        U[i, i] += γt
    end
    U
end

function evaluate_and_jacobian!(u, U, H::FixedPointHomotopy, x, t, c::FixedPointHomotopyCache)
    evaluate_and_jacobian!(u, U, H.F, x, c.F)

    u .= (γ(H) * t) .* (x .- H.x₀) .+ (1 - t) .* u
    U .= (1 - t) .* U
    γt = γ(H) * t
    for i=1:length(x)
        U[i, i] += γt
    end
    U

    nothing
end

function jacobian_and_dt!(U, u, H::FixedPointHomotopy, x, t, c::FixedPointHomotopyCache)
    evaluate_and_jacobian!(u, U, H.F, x, c.F)

    U .= (1 - t) .* U
    γt = γ(H) * t
    for i=1:length(x)
        U[i, i] += γt
    end
    u .= γ(H) .* (x .- H.x₀) .- u

    nothing
end
