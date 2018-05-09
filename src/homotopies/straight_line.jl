export StraightLineHomotopy, StraightLineHomotopyCache

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
(H::StraightLineHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

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

    u = Systems.evaluate(H.start, x, start_cache)
    U = Systems.jacobian(H.start, x, start_cache)

    StraightLineHomotopyCache(start_cache, target_cache, u, U)
end

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
    Systems.evaluate!(c.u, H.start, x, c.start)
    Systems.evaluate!(u, H.target, x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u

    u
end

function dt!(u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate!(c.u, H.start, x, c.start)
    Systems.evaluate!(u, H.target, x, c.target)

    u .= γ(H) .* c.u .- u

    u
end

function jacobian!(U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.jacobian!(c.U, H.start, x, c.start)
    Systems.jacobian!(U, H.target, x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    U
end

function evaluate_and_jacobian!(u, U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, H.start, x, c.start)
    Systems.evaluate_and_jacobian!(u, U, H.target, x, c.target)

    u .= (γ(H) * t) .* c.u .+ (1 - t) .* u
    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U

    nothing
end

function jacobian_and_dt!(U, u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    Systems.evaluate_and_jacobian!(c.u, c.U, H.start, x, c.start)
    Systems.evaluate_and_jacobian!(u, U, H.target, x, c.target)

    U .= (γ(H) * t) .* c.U .+ (1 - t) .* U
    u .= γ(H) .* c.u .- u

    nothing
end

function evaluate(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    M, N = size(H)
    evaluate!(similar(c.u, M), H, x, t, c)
end
function jacobian(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    M, N = size(H)
    jacobian!(similar(c.U, M, N), H, x, t, c)
end
function dt(H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    M, N = size(H)
    dt!(similar(c.u, M), H, x, t, c)
end
