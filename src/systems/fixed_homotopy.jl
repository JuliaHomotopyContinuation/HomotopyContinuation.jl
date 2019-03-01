export FixedHomotopy

"""
    FixedHomotopy(H, t) <: AbstractSystem

Fix a homotopy `H(x,t)` at `t`
"""
mutable struct FixedHomotopy{Hom<:AbstractHomotopy, T} <: AbstractSystem
    H::Hom
    t::T
end

struct FixedHomotopyCache{HC<:AbstractHomotopyCache} <: AbstractSystemCache
    cache::HC
end

cache(FH::FixedHomotopy, x) = FixedHomotopyCache(cache(FH.H, x, FH.t))

Base.size(F::FixedHomotopy) = size(F.H)

evaluate!(u, F::FixedHomotopy, x, c::FixedHomotopyCache) = evaluate!(u, F.H, x, F.t, c.cache)
evaluate(F::FixedHomotopy, x, c::FixedHomotopyCache) = evaluate(F.H, x, F.t, c.cache)
evaluate!(u, F::FixedHomotopy, x) = evaluate!(u, F.H, x, F.t)
evaluate(F::FixedHomotopy, x) = evaluate(F.H, x, F.t)

jacobian!(U, F::FixedHomotopy, x, c::FixedHomotopyCache) = jacobian!(U, F.H, x, F.t, c.cache)
jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache) = jacobian(F.H, x, F.t, c.cache)
jacobian!(U, F::FixedHomotopy, x) = jacobian!(U, F.H, x, F.t)
jacobian(F::FixedHomotopy, x) = jacobian(F.H, x, F.t)

evaluate_and_jacobian!(u, U, F::FixedHomotopy, x, c::FixedHomotopyCache) = evaluate_and_jacobian!(u, U, F.H, x, F.t, c.cache)
evaluate_and_jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache) = evaluate_and_jacobian(F.H, x, F.t, c.cache)
evaluate_and_jacobian!(u, U, F::FixedHomotopy, x) = evaluate_and_jacobian!(u, U, F.H, x, F.t)
evaluate_and_jacobian(F::FixedHomotopy, x) = evaluate_and_jacobian(F.H, x, F.t)
