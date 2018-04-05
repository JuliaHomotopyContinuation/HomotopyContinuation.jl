module Systems

import StaticPolynomials
const SP = StaticPolynomials
import MultivariatePolynomials
const MP = MultivariatePolynomials

export AbstractSystem,
    AbstractSystemCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    evaluate_and_jacobian,
    evaluate_and_jacobian!

export NullCache,
    SPSystem

"""
    AbstractSystem{M, N}

Representing a system of `M` functions in `N` variables.
"""
abstract type AbstractSystem{M, N} end

"""
    AbstractSystemCache

A cache to avoid allocations for the evaluation of an [`AbstractSystem`](@ref).
"""
abstract type AbstractSystemCache end

"""
    NullCache

An empty cache if no cache is necessary.
"""
struct NullCache <: AbstractSystemCache end

"""
    cache(F::AbstractSystem, x)

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`. The default implementation returns a `NullCache`.
"""
cache(F::AbstractSystem, x) = NullCache()

"""
    evaluate!(u, F::AbstractSystem, x [, cache::AbstractSystemCache])

Evaluate the system `F` at `x` and store the result in `u`.
"""
evaluate!(u, F::AbstractSystem, x, ::NullCache) = evaluate!(u, F, x)

"""
    evaluate(F::AbstractSystem, x::AbstractVector [, cache::AbstractSystemCache])

Evaluate the system `F` at `x`.
"""
evaluate(F::AbstractSystem, x, ::NullCache) = evaluate(F, x)


"""
    jacobian!(u, F::AbstractSystem, x [, cache::AbstractSystemCache])

Evaluate the Jacobian of the system `F` at `x` and store the result in `u`.
"""
jacobian!(U, F::AbstractSystem, x, ::NullCache) = jacobian!(U, F, x)

"""
    jacobian(F::AbstractSystem, x [, cache::AbstractSystemCache])

Evaluate the Jacobian of the system `F` at `x`.
"""
jacobian(F::AbstractSystem, x, ::NullCache) = jacobian(F, x)


# Optional
"""
    evaluate_and_jacobian!(u, U, F, x [, cache::AbstractSystemCache])

Evaluate the system `F` and its Jacobian at `x` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, S::AbstractSystem, x::AbstractVector, cache)
    evaluate!(u, S, x, cache)
    jacobian!(U, S, x, cache)
    nothing
end

"""
    evaluate_and_jacobian(F::AbstractSystem, x [, cache::AbstractSystemCache])

Evaluate the system `F` and its Jacobian at `x`.
"""
function evaluate_and_jacobian(S::AbstractSystem, x::AbstractVector, cache)
    u = evaluate(S, x, cache)
    U = jacobian(S, x, cache)
    u, U
end

# Derived
Base.size(S::AbstractSystem{M, N}) where {M, N} = (M, N)


# StaticPolynomialsSystem implementation

"""
    SPSystem(polynomials) <: AbstractSystem

Create a system using the `StaticPolynomials` package.
"""
struct SPSystem{M, N, T, S<:SP.AbstractSystem{T, M, N}} <: AbstractSystem{M, N}
    system::S
end

SPSystem(s::S) where {T,M,N,S<:SP.AbstractSystem{T,M,N}} = SPSystem{M, N, T, S}(s)
SPSystem(polys::Vector{<:MP.AbstractPolynomial}) = SPSystem(SP.system(polys))

evaluate!(u, F::SPSystem, x) = SP.evaluate!(u, F.system, x)
evaluate(F::SPSystem, x) = SP.evaluate(F.system, x)
jacobian!(U, F::SPSystem, x) = SP.jacobian!(U, F.system, x)
jacobian(F::SPSystem, x) = SP.jacobian(F.system, x)
evaluate_and_jacobian!(u, U, F::SPSystem, x) = SP.evaluate_and_jacobian!(u, U, F.system, x)
evaluate_and_jacobian(F::SPSystem, x) = SP.evaluate_and_jacobian(F.system, x)

end
