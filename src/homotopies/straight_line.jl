export StraightLineHomotopy, StraightLineHomotopyCache, gamma

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
    start_cache = cache(H.start, x)
    target_cache = cache(H.target, x)
    u = complex.(float.(evaluate(H.start, x, start_cache)))
    U = complex.(float.(jacobian(H.start, x, start_cache)))
    StraightLineHomotopyCache(start_cache, target_cache, u, U)
end

Base.size(H::StraightLineHomotopy) = size(H.start)

"""
    gamma(H::StraightLineHomotopy)

Obtain the gamma used in the StraightLineHomotopy.
"""
gamma(H::StraightLineHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the StraightLineHomotopy.
"""
γ(H::StraightLineHomotopy) = gamma(H)

function evaluate!(u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    @inbounds evaluate!(c.u, H.start, x, c.start)
    @inbounds evaluate!(u, H.target, x, c.target)

    @inbounds for i in eachindex(u)
        u[i] = (γ(H) * t) * c.u[i] + (1.0 - t) * u[i]
    end

    u
end

function dt!(u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    @inbounds evaluate!(c.u, H.start, x, c.start)
    @inbounds evaluate!(u, H.target, x, c.target)

    @inbounds for i in eachindex(u)
        u[i] = γ(H) * c.u[i] - u[i]
    end
    u
end

function jacobian!(U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    @inbounds jacobian!(c.U, H.start, x, c.start)
    @inbounds jacobian!(U, H.target, x, c.target)

    @inbounds for i in eachindex(U)
        U[i] = (γ(H) * t) * c.U[i] + (1.0 - t) * U[i]
    end

    U
end

function evaluate_and_jacobian!(u, U, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    @inbounds evaluate_and_jacobian!(c.u, c.U, H.start, x, c.start)
    @inbounds evaluate_and_jacobian!(u, U, H.target, x, c.target)

    @inbounds for i in eachindex(u)
        u[i] = (γ(H) * t) * c.u[i] + (1.0 - t) * u[i]
    end
    @inbounds for i in eachindex(U)
        U[i] = (γ(H) * t) * c.U[i] + (1.0 - t) * U[i]
    end

    nothing
end


function jacobian_and_dt!(U, u, H::StraightLineHomotopy, x, t, c::StraightLineHomotopyCache)
    @inbounds evaluate_and_jacobian!(c.u, c.U, H.start, x, c.start)
    @inbounds evaluate_and_jacobian!(u, U, H.target, x, c.target)

    @inbounds for i in eachindex(U)
        U[i] = (γ(H) * t) * c.U[i] + (1.0 - t) * U[i]
    end
    @inbounds for i in eachindex(u)
        u[i] = γ(H) * c.u[i] - u[i]
    end

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
