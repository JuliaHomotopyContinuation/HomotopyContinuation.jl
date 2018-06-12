module SystemsBase

export AbstractSystem,
    AbstractSystemCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    evaluate_and_jacobian,
    evaluate_and_jacobian!

export NullCache

"""
    AbstractSystem

Representing a system of polynomials.
"""
abstract type AbstractSystem end

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
function evaluate_and_jacobian!(u, U, F::AbstractSystem, x, cache)
    evaluate!(u, F, x, cache)
    jacobian!(U, F, x, cache)
    nothing
end

"""
    evaluate_and_jacobian(F::AbstractSystem, x [, cache::AbstractSystemCache])

Evaluate the system `F` and its Jacobian at `x`.
"""
function evaluate_and_jacobian(F::AbstractSystem, x, cache)
    u = evaluate(F, x, cache)
    U = jacobian(F, x, cache)
    u, U
end

Base.size(F::AbstractSystem, i::Integer) = size(F)[i]
Base.length(F::AbstractSystem) = size(F, 1)

end
