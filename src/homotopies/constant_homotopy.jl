export ConstantHomotopy

"""
    ConstantHomotopy(F::AbstractSystem)

Intepret the system `F` as the (constant) homotopy ``H(x,t)=F(x)``.
"""
struct ConstantHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S
end

Base.size(H::ConstantHomotopy) = size(H.system)

struct ConstantHomotopyCache{HC<:AbstractSystemCache} <: AbstractHomotopyCache
    cache::HC
end

cache(H::ConstantHomotopy, x, t) = ConstantHomotopyCache(cache(H.system, x))

const CH = ConstantHomotopy
const CHC = ConstantHomotopyCache

@propagate_inbounds evaluate!(u, H::CH, x, t, c::CHC) = evaluate!(u, H.system, x, c.cache)
@propagate_inbounds jacobian!(U, H::CH, x, t, c::CHC) = jacobian!(U, H.system, x, c.cache)
@propagate_inbounds function dt!(u, H::CH, x, t, c::CHC)
    u .= zero(eltype(u))
    u
end
@propagate_inbounds function evaluate_and_jacobian!(u, U, H::CH, x, t, c::CHC)
    evaluate_and_jacobian!(u, U, H.system, x, c.cache)
end
@propagate_inbounds function jacobian_and_dt!(U, u, H::CH, x, t, c::CHC)
    jacobian!(U, H.system, x, c.cache)
    u .= zero(eltype(u))
    nothing
end
evaluate(H::CH, x, t, c::CHC) = evaluate(H.system, x, c.cache)
jacobian(H::CH, x, t, c::CHC) = jacobian(H.system, x, c.cache)
function dt(H::CH, x, t, c::CHC)
    u = evaluate(H.system, x, c.cache)
    u .= zero(eltype(u))
    u
end
