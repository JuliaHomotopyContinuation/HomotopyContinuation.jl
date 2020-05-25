export TraceTestSystem

"""
    TraceTestSystem(F, p₀, v) <: AbstractSystem

The system `G(x,l) := [F(x, p₀+x[end]v); l ⋅[x[1:end-1];1])`.
"""
struct TraceTestSystem{S<:AbstractSystem,T} <: AbstractSystem
    F::S
    p₀::Vector{T}
    v::Vector{T}
end

function TraceTestSystem(F, p::AbstractVector, v::AbstractVector)
    TraceTestSystem(F, promote(p, v)...)
end

struct TraceTestSystemCache{SC<:AbstractSystemCache,T1,T2} <: AbstractSystemCache
    u::Vector{T1}
    U::Matrix{T1}
    Up::Matrix{T1}
    p::Vector{T2}
    cache::SC
end

function cache(F::TraceTestSystem, x, l)
    c = cache(F.F, x, l)
    u = evaluate(F.F, x, l, c)
    U = similar(u, size(F.F))
    Up = similar(u, size(F.F)[1], length(F.p₀))
    p = F.p₀ + 0.1 * F.v
    TraceTestSystemCache(u, U, Up, p, c)
end

Base.size(F::TraceTestSystem) = size(F.F) .+ 1

function evaluate!(u, F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    λ = x[end]
    for i in eachindex(c.p)
        c.p[i] = F.p₀[i] + λ * F.v[i]
    end
    evaluate!(c.u, F.F, x, c.p, c.cache)
    n = length(x)
    for i = 1:length(c.u)
        u[i] = c.u[i]
    end
    u[end] = l[n]
    for i = 1:n-1
        u[end] += l[i] * x[i]
    end

    u
end

function jacobian!(U, F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    λ = x[end]
    for i in eachindex(c.p)
        c.p[i] = F.p₀[i] + λ * F.v[i]
    end
    jacobian!(c.U, F.F, x, c.p, c.cache)
    differentiate_parameters!(c.Up, F.F, x, c.p, c.cache)
    n, m = length(x), size(U, 1)
    for j = 1:n-1, i = 1:m-1
        U[i, j] = c.U[i, j]
    end
    for i = 1:m
        U[i, n] = zero(eltype(U))
    end
    for k = 1:size(c.Up, 2), i = 1:m-1
        U[i, n] += c.Up[i, k] * F.v[k]
    end

    for j = 1:n-1
        U[m, j] = l[j]
    end
    U
end


function evaluate_and_jacobian!(u, U, F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    λ = x[end]
    for i in eachindex(c.p)
        c.p[i] = F.p₀[i] + λ * F.v[i]
    end
    n, m = length(x), size(U, 1)

    evaluate_and_jacobian!(c.u, c.U, F.F, x, c.p, c.cache)
    for i = 1:length(c.u)
        u[i] = c.u[i]
    end
    u[end] = l[n]
    for i = 1:n-1
        u[end] += l[i] * x[i]
    end

    differentiate_parameters!(c.Up, F.F, x, c.p, c.cache)

    for j = 1:n-1, i = 1:m-1
        U[i, j] = c.U[i, j]
    end
    for i = 1:m
        U[i, n] = zero(eltype(U))
    end
    for k = 1:size(c.Up, 2), i = 1:m-1
        U[i, n] += c.Up[i, k] * F.v[k]
    end

    for j = 1:n-1
        U[m, j] = l[j]
    end

    nothing
end

function differentiate_parameters!(U, F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    U .= zero(eltype(U))
    for i = 1:length(l)-1
        U[end, i] = x[i]
    end
    U[end, length(l)] = one(eltype(U))
    U
end

function differentiate_parameters(F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    differentiate_parameters!(similar(c.U, size(F)[1], length(l)), F, x, l, c)
end

function evaluate(F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    evaluate!(similar(c.u, size(F)[1]), F, x, l, c)
end
function jacobian(F::TraceTestSystem, x, l, c::TraceTestSystemCache)
    jacobian!(similar(c.U, size(F)), F, x, l, c)
end
