module Correctors

using ..Homotopies
using ..Utilities

export AbstractCorrector,
    AbstractCorrectorCache,
    Result,
    cache,
    correct!,
    Newton, NewtonCache

# Interface
abstract type AbstractCorrector end
abstract type AbstractCorrectorCache end

"""
    Result{T}

Structure holding information about a `correct!` step. The fields are
* `converged::Bool` Indicating whether the correction was successfull
* `res::T` The residual of the final result.
* `iters::Int` The number of iterations used.
"""
struct Result{T}
    converged::Bool
    res::T
    iters::Int
end

"""
    cache(::AbstractCorrector, ::HomotopyWithCache{M, N}, x, t)::AbstractCorrectorCache

Construct a cache to avoid allocations.
"""
function cache end


"""
    correct!(xnext, ::AbstractCorrector, ::AbstractCorrectorCache, H::HomotopyWithCache, x, t, tol)::Result

Perform a correction step such that in the end `H(xnext,t) < tol`.
Returns a [`Result`](@ref).
"""
function correct! end


# Newton
"""
    Newton()

A classical simple Newton operator for square linear systems using the LU factorization
to solve the linear systems.
"""
struct Newton <: AbstractCorrector end

struct NewtonCache{T} <: AbstractCorrectorCache
    A::Matrix{T}
    b::Vector{T}
end

function cache(::Newton, H::HomotopyWithCache, x, t)
    A = jacobian(H, x, t)
    b = evaluate(H, x, t)
    NewtonCache(A, b)
end

function correct!(xnext, ::Newton, cache::NewtonCache, H::HomotopyWithCache, x, t, tol, maxiters)
    A, b = cache.A, cache.b
    evaluate_and_jacobian!(b, A, H, x, t)
    solve_with_lu_inplace!(A, b)
    @. xnext = xnext - b

    k = 1
    while true
        evaluate!(b, H, xnext, t)

        res = infinity_norm(b)
        if res < tol
            return Result(true, res, k)
        elseif k ≥ maxiters
            return Result(false, res, k)
        end

        k += 1

        # put jacobian in A
        jacobian!(A, H, xnext, t)

        solve_with_lu_inplace!(A, b)
        @. xnext = xnext - b
    end
end

end
