module Homotopies

import ..HomotopiesBase: AbstractHomotopy,
    AbstractHomotopyCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    precondition!, update!

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import ..Utilities: randomish_gamma

export HomotopyWithCache,
    AbstractHomotopy,
    AbstractHomotopyCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    precondition!, update!

struct HomotopyWithCache{H<:AbstractHomotopy, C<:AbstractHomotopyCache} <: AbstractHomotopy
    homotopy::H
    cache::C
end

HomotopyWithCache(H::AbstractHomotopy, x, t) = HomotopyWithCache(H, cache(H, x, t))

(H::HomotopyWithCache)(x, t) = evaluate(H, x, t)

evaluate(H::HomotopyWithCache, x, t) = evaluate(H.homotopy, x, t, H.cache)
evaluate!(u, H::HomotopyWithCache, x, t) = evaluate!(u, H.homotopy, x, t, H.cache)
jacobian(H::HomotopyWithCache, x, t) = jacobian(H.homotopy, x, t, H.cache)
jacobian!(U, H::HomotopyWithCache, x, t) = jacobian!(U, H.homotopy, x, t, H.cache)
dt(H::HomotopyWithCache, x, t) = dt(H.homotopy, x, t, H.cache)
dt!(u, H::HomotopyWithCache, x, t) = dt!(u, H.homotopy, x, t, H.cache)
jacobian_and_dt(H::HomotopyWithCache, x, t) = jacobian_and_dt(H.homotopy, x, t, H.cache)
jacobian_and_dt!(U, u, H::HomotopyWithCache, x, t) = jacobian_and_dt!(U, u, H.homotopy, x, t, H.cache)
evaluate_and_jacobian(H::HomotopyWithCache, x, t) = evaluate_and_jacobian(H.homotopy, x, t, H.cache)
evaluate_and_jacobian!(u, U, H::HomotopyWithCache, x, t) = evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.cache)
precondition!(H, x, t) = precondition!(H.homotopy, x, t, H.cache)
update!(H, x, t) = update!(H.homotopy, x, t, H.cache)

include("homotopies/straight_line.jl")
include("homotopies/fixed_point.jl")
include("homotopies/patched_homotopy.jl")
end
