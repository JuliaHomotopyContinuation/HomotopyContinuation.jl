export AbstractSystem,
    AbstractSystemCache,
    SystemNullCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    differentiate_parameters!,
    differentiate_parameters,
    evaluate_and_jacobian,
    evaluate_and_jacobian!


"""
    SystemNullCache

An empty cache if no cache is necessary.
"""
struct SystemNullCache <: AbstractSystemCache end

"""
    cache(F::AbstractSystem, x)::AbstractSystemCache

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`.

    cache(F::AbstractSystem, x, p)::AbstractSystemCache

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x` and parameters `p`.
"""
cache(F::AbstractSystem, args...) = throw(MethodError(cache, tuple(F, args...)))


"""
    evaluate!(u, F::AbstractSystem, x, cache::AbstractSystemCache)

Evaluate the system `F` at `x` and store the result in `u`.

    evaluate!(u, F::AbstractSystem, x, p, cache::AbstractSystemCache)

Evaluate the system `F` at `x` and parameters `p` and store the result in `u`.
"""
evaluate!(u, F::AbstractSystem, args...) = error(MethodError(evaluate!, tuple(u, F, args...)))

"""
    evaluate(F::AbstractSystem, x::AbstractVector, cache=cache(F, x))

Evaluate the system `F` at `x`.

    evaluate(F::AbstractSystem, x::AbstractVector, p, cache=cache(F, x))

Evaluate the system `F` at `x` and parameters `p`.
"""
evaluate(F::AbstractSystem, x, c::AbstractSystemCache=cache(F, x)) = evaluate(F, x, c)
evaluate(F::AbstractSystem, x, p, c::AbstractSystemCache=cache(F,x,p)) = evaluate(F, x, p, c)


"""
    jacobian!(u, F::AbstractSystem, x , cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and store the result in `u`.

    jacobian!(u, F::AbstractSystem, x , p, cache::AbstractSystemCache)

Evaluate the Jacobian of the system `F` at `x` and parameters `p` and store the result in `u`.
"""
jacobian!(u, F::AbstractSystem, args...) = error(MethodError(jacobian!, tuple(u, F, args...)))

"""
    jacobian(F::AbstractSystem, x, cache=cache(F, x))

Evaluate the Jacobian of the system `F` at `x`.

    jacobian(F::AbstractSystem, x , p, cache=cache(F, x))

Evaluate the Jacobian of the system `F` at `x` and parameters `p`.
"""
jacobian(F::AbstractSystem, x, c::AbstractSystemCache=cache(F, x)) = jacobian(F, x, c)
jacobian(F::AbstractSystem, x, p, c::AbstractSystemCache=cache(F,x,p)) = jacobian(F, x, p, c)


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
differentiate_parameters(F::AbstractSystem, x, p, c::AbstractSystemCache=cache(F, x)) = differentiate_parameters(F, x, p, c)

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
function evaluate_and_jacobian(F::AbstractSystem, x, c::AbstractSystemCache=cache(F, x))
    u = evaluate(F, x, c)
    U = jacobian(F, x, c)
    u, U
end
"""
    evaluate_and_jacobian(F::AbstractSystem, x, p, cache=cache(F, x))

Evaluate the system `F` and its Jacobian at `x` and parameters `p`.
"""
function evaluate_and_jacobian(F::AbstractSystem, x, p, c::AbstractSystemCache=cache(F, x))
    u = evaluate(F, x, p, c)
    U = jacobian(F, x, p, c)
    u, U
end

"""
    Base.size(F::AbstractSystem)

Returns a tuple `(m, n)` indicating that `F` is a system of `m` polynomials `m` in `n` variables.
"""
Base.size(::AbstractSystem) = error("Mandatory to define `Base.size` for `AbstractSystem`s")
Base.size(F::AbstractSystem, i::Integer) = size(F)[i]
Base.length(F::AbstractSystem) = size(F, 1)

include("systems/sp_system.jl")
include("systems/fixed_homotopy.jl")
include("systems/fp_system.jl")
include("systems/totaldegree_system.jl")
include("systems/multihom_totaldegree_system.jl")
include("systems/composition_system.jl")
include("systems/fixed_parameter_system.jl")
include("systems/patched_system.jl")
include("systems/squared_up_system.jl")
