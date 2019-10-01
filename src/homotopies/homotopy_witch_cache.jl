export HomotopyWithCache

struct HomotopyWithCache{H<:AbstractHomotopy,C<:AbstractHomotopyCache} <: AbstractHomotopy
    homotopy::H
    cache::C
end

HomotopyWithCache(H::AbstractHomotopy, x, t) = HomotopyWithCache(H, cache(H, x, t))

cache(::HomotopyWithCache, x, t) = HomotopyNullCache()

(H::HomotopyWithCache)(x, t) = evaluate(H, x, t)

Base.size(H::HomotopyWithCache) = size(H.homotopy)
evaluate(H::HomotopyWithCache, x, t, ::HNC = HNC()) = evaluate(H.homotopy, x, t, H.cache)
@propagate_inbounds evaluate!(u, H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    evaluate!(u, H.homotopy, x, t, H.cache)
jacobian(H::HomotopyWithCache, x, t, ::HNC = HNC()) = jacobian(H.homotopy, x, t, H.cache)
@propagate_inbounds jacobian!(U, H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    jacobian!(U, H.homotopy, x, t, H.cache)
dt(H::HomotopyWithCache, x, t, ::HNC = HNC()) = dt(H.homotopy, x, t, H.cache)
@propagate_inbounds dt!(u, H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    dt!(u, H.homotopy, x, t, H.cache)
jacobian_and_dt(H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    jacobian_and_dt(H.homotopy, x, t, H.cache)
@propagate_inbounds jacobian_and_dt!(U, u, H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    jacobian_and_dt!(U, u, H.homotopy, x, t, H.cache)
evaluate_and_jacobian(H::HomotopyWithCache, x, t, ::HNC = HNC()) =
    evaluate_and_jacobian(H.homotopy, x, t, H.cache)
@propagate_inbounds evaluate_and_jacobian!(
    u,
    U,
    H::HomotopyWithCache,
    x,
    t,
    ::HNC = HNC(),
) = evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.cache)
basehomotopy(H::HomotopyWithCache) = basehomotopy(H.homotopy)
