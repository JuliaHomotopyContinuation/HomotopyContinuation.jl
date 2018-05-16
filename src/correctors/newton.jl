export Newton

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
    ldiv_lu!(A, b)
    @. xnext = x - b
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

        ldiv_lu!(A, b)
        @. xnext = xnext - b
    end
end
