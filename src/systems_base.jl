module SystemsBase

export AbstractSystem,
    AbstractSystemCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    differentiate_parameters!,
    differentiate_parameters,
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
    cache(F::AbstractSystem, x)::AbstractSystemCache

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`.

    cache(F::AbstractSystem, x, p)::AbstractSystemCache

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x` and parameters `p`.
"""
function cache end


"""
    evaluate!(u, F::AbstractSystem, x, cache::AbstractSystemCache)

Evaluate the system `F` at `x` and store the result in `u`.

    evaluate!(u, F::AbstractSystem, x, p, cache::AbstractSystemCache)

Evaluate the system `F` at `x` and parameters `p` and store the result in `u`.
"""
function evaluate! end

"""
    evaluate(F::AbstractSystem, x::AbstractVector, cache=cache(F, x))

Evaluate the system `F` at `x`.

    evaluate(F::AbstractSystem, x::AbstractVector, p, cache=cache(F, x))

Evaluate the system `F` at `x` and parameters `p`.
"""
evaluate(F::AbstractSystem, x, c=cache(F, x)) = evaluate(F, x, c)


"""
    jacobian!(u, F::AbstractSystem, x , cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and store the result in `u`.

    jacobian!(u, F::AbstractSystem, x , p, cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and parameters `p` and store the result in `u`.
"""
function jacobian! end

"""
    jacobian(F::AbstractSystem, x, cache=cache(F, x))

Evaluate the Jacobian of the system `F` at `x`.

    jacobian(F::AbstractSystem, x , p, cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and parameters `p`.
"""
jacobian(F::AbstractSystem, x, c=cache(F, x)) = jacobian(F, x, c)


"""
    differentiate_parameters!(u, F::AbstractSystem, x, p, cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and parameters `p` w.r.t. the parameters
and store the result in `u`.
"""
function differentiate_parameters! end

"""
    differentiate_parameters(F::AbstractSystem, x, p, cache=cache(F, x))

Evaluate the Jacobian of the system `F` at `x` and parameters `p` w.r.t. the parameters
"""
differentiate_parameters(F::AbstractSystem, x, c=cache(F, x)) = differentiate_parameters(F, x, c)

# Optional
"""
    evaluate_and_jacobian!(u, U, F, x , cache::AbstractSystemCache)

Evaluate the system `F` and its Jacobian at `x` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, F::AbstractSystem, x, cache::AbstractSystemCache)
    evaluate!(u, F, x, cache)
    jacobian!(U, F, x, cache)
    nothing
end

"""
    evaluate_and_jacobian!(u, U, F, x, p, cache::AbstractSystemCache)

Evaluate the system `F` and its Jacobian at `x` and parameters `p` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, F::AbstractSystem, x, p, cache::AbstractSystemCache)
    evaluate!(u, F, x, p, cache)
    jacobian!(U, F, x, p, cache)
    nothing
end


"""
    evaluate_and_jacobian(F::AbstractSystem, x, cache=cache(F, x))

Evaluate the system `F` and its Jacobian at `x`.
"""
function evaluate_and_jacobian(F::AbstractSystem, x, c=cache(F, x))
    u = evaluate(F, x, c)
    U = jacobian(F, x, c)
    u, U
end
"""
    evaluate_and_jacobian(F::AbstractSystem, x, p, cache=cache(F, x))

Evaluate the system `F` and its Jacobian at `x` and parameters `p`.
"""
function evaluate_and_jacobian(F::AbstractSystem, x, p, c=cache(F, x))
    u = evaluate(F, x, p, c)
    U = jacobian(F, x, p, c)
    u, U
end

Base.size(F::AbstractSystem, i::Integer) = size(F)[i]
Base.length(F::AbstractSystem) = size(F, 1)

end
