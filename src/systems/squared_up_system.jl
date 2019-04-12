export SquaredUpSystem

"""
    SquaredUpSystem(F, A) <: AbstractSystem

Square up a system `F(x)`. Assume `m, n = size(F)`. Then this computes with with the system
``(f_1,...,f_n) + A(f_m+1,...,f_n)``  where ``A`` is a random ``n × (m - n)`` matrix.
"""
struct SquaredUpSystem{S<:AbstractSystem, T} <: AbstractSystem
    F::S
    A::Matrix{T}
end

struct SquaredUpSystemCache{SC<:AbstractSystemCache, T} <: AbstractSystemCache
    u::Vector{T}
    U::Matrix{T}
    cache::SC
end

function cache(F::SquaredUpSystem, x)
    c = cache(F.F, x)
    u, U = evaluate_and_jacobian(F.F, x, c)
    SquaredUpSystemCache(u, U, c)
end

Base.size(F::SquaredUpSystem) = (size(F.A, 1), size(F.F)[2])

function evaluate!(u, F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    evaluate!(c.u, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    for i in 1:n
        u[i] = c.u[i]
        for j in 1:(m-n)
            u[i] += F.A[i, j] * c.u[n + j]
        end
    end
    u
end

function jacobian!(U, F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    jacobian!(c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    n′ = size(F.F)[2]
    for j in 1:n′, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j]
        end
    end
    U
end

function evaluate_and_jacobian!(u, U, F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    evaluate_and_jacobian!(c.u, c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    n′ = size(F.F)[2]

    for i in 1:n
        u[i] = c.u[i]
        for j in 1:(m-n)
            u[i] += F.A[i, j] * c.u[n + j]
        end
    end
    for j in 1:n′, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j]
        end
    end
    nothing
end

function evaluate(F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    u = Vector{_return_type(F,c)}(undef, size(F.A, 1))
    evaluate!(u, F, x, c)
end
function jacobian(F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    U = Matrix{_return_type(F,c)}(undef, size(F.A, 1), size(F.F)[2])
    jacobian!(U, F, x, c)
end

function evaluate_and_jacobian(F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    T = _return_type(F,c)
    u = Vector{T}(undef, size(F.A, 1))
    U = Matrix{T}(undef, size(F.A, 1), size(F.F)[2])
    evaluate_and_jacobian!(u, U, F, x, c)
    u, U
end

_return_type(F::SquaredUpSystem, c::SquaredUpSystemCache) = typeof(F.A[1,1] * c.u[1])
