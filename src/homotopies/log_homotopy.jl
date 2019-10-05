export LogHomotopy

"""
    LogHomotopy(H::AbstractHomotopy)

Reparameterizes the homotopy ``H(x,t)`` as ``Ĥ(x,s) := H(x,e^{-s})``.
"""
struct LogHomotopy{H<:AbstractHomotopy} <: AbstractHomotopy
    homotopy::H
end

Base.size(H::LogHomotopy) = size(H.homotopy)
basehomotopy(H::LogHomotopy) = basehomotopy(H.homotopy)

struct LogHomotopyCache{HC} <: AbstractHomotopyCache
    cache::HC
end

function cache(H::LogHomotopy, x, s)
    t = exp(-s)
    c = cache(H.homotopy, x, t)
    LogHomotopyCache(c)
end

@propagate_inbounds function evaluate!(u, H::LogHomotopy, x, s, c::LogHomotopyCache)
    evaluate!(u, H.homotopy, x, exp(-s), c.cache)
end
@propagate_inbounds function jacobian!(U, H::LogHomotopy, x, s, c::LogHomotopyCache)
    jacobian!(U, H.homotopy, x, exp(-s), c.cache)
end
@propagate_inbounds function dt!(u, H::LogHomotopy, x, s, c::LogHomotopyCache)
    t = exp(-s)
    dt!(u, H.homotopy, x, t, c.cache)
    LinearAlgebra.rmul!(u, -t)
    u
end

@propagate_inbounds function evaluate_and_jacobian!(u, U, H::LogHomotopy, x, s, c::LogHomotopyCache)
    evaluate_and_jacobian!(u, U, H.homotopy, x, exp(-s), c.cache)
end

@propagate_inbounds function jacobian_and_dt!(U, u, H::LogHomotopy, x, s, c::LogHomotopyCache)
    t = exp(-s)
    jacobian_and_dt!(U, u, H.homotopy, x, t, c.cache)
    LinearAlgebra.rmul!(u, -t)

    nothing
end

function evaluate(H::LogHomotopy, x, s, c::LogHomotopyCache)
    evaluate(H.homotopy, x, exp(-s), c.cache)
end
function jacobian(H::LogHomotopy, x, s, c::LogHomotopyCache)
    jacobian(H.homotopy, x, exp(-s), c.cache)
end
function dt(H::LogHomotopy, x, s, c::LogHomotopyCache)
    t = exp(-s)
    u = dt(H.homotopy, x, t, c.cache)
    LinearAlgebra.rmul!(u, -t)
    u
end
