export FixedParameterSystem

"""
    FixedParameterSystem(F, p) <: AbstractSystem

Fix a parametrized system `F(x; p)` at `p`, i.e., it is treated as a system without parameters.
"""
mutable struct FixedParameterSystem{S<:AbstractSystem,P<:AbstractVector} <: AbstractSystem
    F::S
    p::P
end

struct FixedParameterSystemCache{SC<:AbstractSystemCache} <: AbstractSystemCache
    cache::SC
end

cache(F::FixedParameterSystem, x) = FixedParameterSystemCache(cache(F.F, x, F.p))

Base.size(F::FixedParameterSystem) = size(F.F)

evaluate!(u, F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    evaluate!(u, F.F, x, F.p, c.cache)
evaluate(F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    evaluate(F.F, x, F.p, c.cache)
jacobian!(U, F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    jacobian!(U, F.F, x, F.p, c.cache)
jacobian(F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    jacobian(F.F, x, F.p, c.cache)
evaluate_and_jacobian!(u, U, F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    evaluate_and_jacobian!(u, U, F.F, x, F.p, c.cache)
evaluate_and_jacobian(F::FixedParameterSystem, x, c::FixedParameterSystemCache) =
    evaluate_and_jacobian(F.F, x, F.p, c.cache)

"""
    set_parameters!(F::FixedParameterSystem, p)

Update the parameters `p` of `F`.
"""
set_parameters!(F::FixedParameterSystem, p) = begin
    F.p .= p
    F
end
