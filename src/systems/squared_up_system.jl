export SquaredUpSystem

"""
    SquaredUpSystem(F, A, degrees) <: AbstractSystem

Square up a system `F(x)`. Assume `m, n = size(F)`. Then this computes with with the system
``(f_1,...,f_n) + A(f_m+1,...,f_n)``  where ``A`` is a random ``n × (m - n)`` matrix.
"""
struct SquaredUpSystem{S<:AbstractSystem, T} <: AbstractSystem
    F::S
    A::Matrix{T}
    degrees::Vector{Int}
end

function SquaredUpSystem(F::AbstractSystem, A::Matrix)
    SquaredUpSystem(F, A, zeros(Int, sum(size(A))))
end

struct SquaredUpSystemCache{SC<:AbstractSystemCache, T} <: AbstractSystemCache
    u::Vector{T}
    U::Matrix{T}
    cache::SC
    degree_diffs::Matrix{Int}
end

function cache(F::SquaredUpSystem, x)
    c = cache(F.F, x)
    u, U = evaluate_and_jacobian(F.F, x, c)
    n = size(F.A, 1)
    degree_diffs = [F.degrees[i] - F.degrees[n+j] for i=1:n, j=1:size(F.A, 2)]
    SquaredUpSystemCache(u, U, c, degree_diffs)
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

function evaluate!(u, F::SquaredUpSystem, x::PVector{<:Number,1}, c::SquaredUpSystemCache)
    evaluate!(c.u, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    for i in 1:n
        u[i] = c.u[i]
    end
    for j in 1:(m-n), i in 1:n
        u[i] += F.A[i, j] * c.u[n + j] * x[end]^c.degree_diffs[i,j]
    end
    u
end


function jacobian!(U, F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    jacobian!(c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    for j in 1:n, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j]
        end
    end
    U
end

function jacobian!(U, F::SquaredUpSystem, x::PVector{<:Number,1}, c::SquaredUpSystemCache)
    evaluate_and_jacobian!(c.u, c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]
    for j in 1:n+1, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j] * x[end]^c.degree_diffs[i,k]
        end
    end

    # we have to take care of the product rule in the derivative wrt x[n+1]
    for i in 1:n, k in 1:(m - n)
        if c.degree_diffs[i, k] == 1
            U[i, n+1] += F.A[i, k] * c.u[n+k]
        elseif c.degree_diffs[i, k] > 1
            U[i, n+1] += c.degree_diffs[i, k] * F.A[i, k] * c.u[n+k] * x[n+1]^(c.degree_diffs[i, k] - 1)
        end
    end
    U
end

function evaluate_and_jacobian!(u, U, F::SquaredUpSystem, x, c::SquaredUpSystemCache)
    evaluate_and_jacobian!(c.u, c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]

    for i in 1:n
        u[i] = c.u[i]
        for j in 1:(m-n)
            u[i] += F.A[i, j] * c.u[n + j]
        end
    end
    for j in 1:n, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j]
        end
    end
    nothing
end

function evaluate_and_jacobian!(u, U, F::SquaredUpSystem, x::PVector{<:Number, 1}, c::SquaredUpSystemCache)
    evaluate_and_jacobian!(c.u, c.U, F.F, x, c.cache)
    n, m = size(F.A, 1), size(F.F)[1]

    for i in 1:n
        u[i] = c.u[i]
    end
    for j in 1:(m-n), i in 1:n
        u[i] += F.A[i, j] * c.u[n + j] * x[end]^c.degree_diffs[i,j]
    end

    for j in 1:n+1, i in 1:n
        U[i, j] = c.U[i, j]
        # + A * F̂
        for k in 1:(m - n)
            U[i, j] += F.A[i, k] * c.U[n+k, j] * x[end]^c.degree_diffs[i,k]
        end
    end

    # we have to take care of the product rule in the derivative wrt x[n+1]
    for i in 1:n, k in 1:(m - n)
        if c.degree_diffs[i,k] == 1
            U[i, n+1] += F.A[i, k] * c.u[n+k]
        elseif c.degree_diffs[i,k] > 1
            U[i, n+1] += c.degree_diffs[i, k] * F.A[i, k] * c.u[n+k] * x[n+1]^(c.degree_diffs[i,k] - 1)
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
