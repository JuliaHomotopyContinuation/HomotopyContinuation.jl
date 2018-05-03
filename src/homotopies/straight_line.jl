"""
    StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct StraightLineHomotopy{S<:AbstractSystem,T<:AbstractSystem} <: AbstractHomotopy
    start::S
    target::T
    gamma::Complex{Float64}

end
function StraightLineHomotopy(start::AbstractSystem, target::AbstractSystem; gamma=randomish_gamma())
    StraightLineHomotopy(start, target, gamma)
end

"""
    StraightLineHomotopyCache

An simple cache for `StartTargetHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct StraightLineHomotopyCache{SC, TC, T1, T2} <: AbstractHomotopyCache
    start::SC
    target::TC

    u::Vector{T1} # for evaluation
    U::Matrix{T2} # for Jacobian
end

function cache(H::StraightLineHomotopy, x, t)
    start_cache = Systems.cache(H.start, x)
    target_cache = Systems.cache(H.target, x)

    u = Systems.evaluate(H.start, x)
    U = Systems.jacobian(H.start, x)

    StraightLineHomotopyCache(start_cache, target_cache, u, U)
end

start(H::StraightLineHomotopy) = H.start
target(H::StraightLineHomotopy) = H.target

Base.size(H::StraightLineHomotopy) = size(H.start)

"""
    gamma(H::StraightLineHomotopy)

Obtain the gamma used in the StraightLineHomotopy.
"""
Base.Math.gamma(H::StraightLineHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the StraightLineHomotopy.
"""
γ(H::StraightLineHomotopy) = gamma(H)

function evaluate!(u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u

    u
end
function evaluate(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    (γ(H) * t) * G + (1 - t) * F
end
(H::StraightLineHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function dt!(u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate!(c.u, start(H), x, c.start)
    Systems.evaluate!(u, target(H), x, c.target)

    u .= γ(H) .* c.u .- u

    u
end
function dt(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    G = Systems.evaluate(start(H), x, c.start)
    F = Systems.evaluate(target(H), x, c.target)
    γ(H) .* G .- F
end

function jacobian!(U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.jacobian!(c.U, start(H), x, c.start)
    Systems.jacobian!(U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    U
end
function jacobian(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    G = Systems.jacobian(start(H), x, c.start)
    F = Systems.jacobian(target(H), x, c.target)
    (γ(H) * t) .* G .+ (1 - t) .* F
end

function evaluate_and_jacobian!(u, U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u
    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    nothing
end

function jacobian_and_dt!(U, u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, start(H), x, c.start)
    Systems.evaluate_and_jacobian!(u, U, target(H), x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U
    u .= γ(H) .* c.u .- u

    nothing
end
