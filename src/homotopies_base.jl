module HomotopiesBase

import ..SystemsBase: AbstractSystem, AbstractSystemCache

export AbstractHomotopy,
    AbstractHomotopyCache,
    NullCache,
    nvariables,
    cache,
    evaluate!,
    evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    precondition!, update!


"""
    AbstractHomotopy

Representing a homotopy.
"""
abstract type AbstractHomotopy end

abstract type AbstractParameterHomotopy <: AbstractHomotopy end


# Cache
"""
    AbstractHomotopyCache

A cache to avoid allocations for the evaluation of an [`AbstractHomotopy`](@ref).
"""
abstract type AbstractHomotopyCache end

"""
    NullCache

The default `AbstractHomotopyCache` containing nothing.
"""
struct NullCache <: AbstractHomotopyCache end

"""
    cache(H::AbstractHomotopy, x, t)

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`.
"""
function cache end
cache(H, x, t) = NullCache()





# Homotopy API

"""
    evaluate!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
function evaluate! end

"""
    evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)`.
"""
evaluate(H::AbstractHomotopy, x, t) = evaluate(H, x, t, cache(H, x, t))
function evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    u = fill(zero(eltype(x)), size(H)[1])
    evaluate!(u, H, x, t, cache)
    u
end

"""
    dt!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
function dt! end

"""
    dt(H::AbstractHomotopy, x::AbstractVector, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)`.
"""
dt(H::AbstractHomotopy, x, t) = dt(H, x, t, cache(H, x, t))
function dt(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    u = fill(zero(eltype(x)), size(H)[1])
    dt!(u, H, x, t, cache)
    u
end

"""
    jacobian!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the Jacobian of the homotopy `H` at `(x, t)` and store the result in `u`.
"""
function jacobian! end

"""
    jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the Jacobian of the homotopy `H` at `(x, t)`.
"""
jacobian(H::AbstractHomotopy, x, t) = jacobian(H, x, t, cache(H, x, t))
function jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    U = fill(zero(eltype(x)), size(H))
    jacobian!(U, H, x, t, cache)
    U
end

# Optional
"""
    evaluate_and_jacobian!(u, U, F, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its Jacobian at `(x, t)` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    evaluate!(u, H, x, t, c)
    jacobian!(U, H, x, t, c)
    nothing
end

"""
    evaluate_and_jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its Jacobian at `(x, t)`.
"""
function evaluate_and_jacobian(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    u = evaluate(H, x, t, c)
    U = jacobian(H, x, t, c)
    u, U
end

"""
    jacobian_and_dt!(U, u, H, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)` and store the results in `u` (evalution)
and `v` (∂t).
"""
function jacobian_and_dt!(U, u, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    jacobian!(U, H, x, t, c)
    dt!(u, H, x, t, c)
    nothing
end

"""
    evaluate_and_dt(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)`.
"""
function jacobian_and_dt(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    U = jacobian(H, x, t, c)
    u = dt(H, x, t, c)
    U, u
end

"""
    nvariables(H::AbstractHomotopy)

Returns the number of variables of the homotopy `H`.
"""
nvariables(H::AbstractHomotopy) = last(size(H))

"""
    precondition!(H::AbstractHomotopy, x, t, cache)

Prepare a homotopy for things like pathtracking starting at `x` and `t`.
This can modify `x` as well as `H` and anything in `cache`.
By default this is a no-op. If `H` wraps another homotopy this should call
`precondition!` on this as well.
"""
precondition!(H::AbstractHomotopy, x, t, cache) = nothing
"""
    update!(H::AbstractHomotopy, x, t, cache)

Update a homotopy for new values of `x` and `x`, i.e., update an affine patch.
This can modify `x` as well as `H` and anything in `cache`.
By default this is a no-op. If `H` wraps another homotopy this should call
`update!` on this as well.
"""
update!(H::AbstractHomotopy, x, t, cache) = nothing

end
