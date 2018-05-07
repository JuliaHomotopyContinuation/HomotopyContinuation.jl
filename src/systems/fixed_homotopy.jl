export FixedHomotopy

const HB = HomotopiesBase

"""
    FixedHomotopy(H, t) <: AbstractSystem

Fix a homotopy `H(x,t)` at `t`
"""
struct FixedHomotopy{Hom<:HB.AbstractHomotopy, T} <: AbstractSystem
    H::Hom
    t::T
end

struct FixedHomotopyCache{HC<:HB.AbstractHomotopyCache} <: AbstractSystemCache
    cache::HC
end

cache(FH::FixedHomotopy, x) = FixedHomotopyCache(HB.cache(FH.H, x, FH.t))

Base.size(F::FixedHomotopy) = size(F.H)

function evaluate!(u, F::FixedHomotopy, x, c::FixedHomotopyCache)
    HB.evaluate!(u, F.H, x, F.t, c.cache)
end
evaluate(F::FixedHomotopy, x, c::FixedHomotopyCache) = HB.evaluate(F.H, x, F.t, c.cache)
function jacobian!(U, F::FixedHomotopy, x, c::FixedHomotopyCache)
    HB.jacobian!(U, F.H, x, F.t, c.cache)
end
function jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache)
    HB.jacobian(F.H, x, F.t, c.cache)
end
function evaluate_and_jacobian!(u, U, F::FixedHomotopy, x, c::FixedHomotopyCache)
    HB.evaluate_and_jacobian!(u, U, F.H, x, F.t, c.cache)
end
function evaluate_and_jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache)
    HB.evaluate_and_jacobian(F.H, x, F.t, c.cache)
end
