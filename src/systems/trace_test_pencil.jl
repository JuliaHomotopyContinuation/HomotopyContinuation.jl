export TraceTestPencil

"""
    TraceTestPencil(F, p₀, v) <: AbstractSystem

The system `H(x,t) := [F(x, p₀+x[end]v); l ⋅[x[1:end-1];1] * λ + τ*γ)`.
"""
struct TraceTestPencil{S<:AbstractSystem,T} <: AbstractSystem
    F::TraceTestSystem{S,T}
    l::Vector{T}
    γ::T
end

function TraceTestPencil(F::TraceTestSystem, l)
    TraceTestPencil(F, l, cis(2π * rand()))
end

struct TraceTestPencilCache{SC<:TraceTestSystemCache} <: AbstractSystemCache
    cache::SC
end

function cache(F::TraceTestPencil, x, τ)
    c = cache(F.F, x, F.l)
    TraceTestPencilCache(c)
end

Base.size(F::TraceTestPencil) = size(F.F)

function evaluate!(u, F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    λ = x[end]
    evaluate!(u, F.F, x, F.l, c.cache)
    u[end] = u[end] * λ + τ[1] * F.γ
    u
end
function evaluate(F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    evaluate!(evaluate(F.F, x, F.l, c.cache), F, x, τ, c)
end

function jacobian!(U, F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    λ = x[end]
    jacobian!(U, F.F, x, F.l, c.cache)
    for j = 1:length(x)-1
        U[end, j] *= λ
    end
    U[end, end] = F.l[end]
    for i = 1:length(x)-1
        U[end, end] += F.l[i] * x[i]
    end
    U
end
function jacobian(F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    jacobian!(jacobian(F.F, x, F.l, c.cache), F, x, τ, c)
end

function evaluate_and_jacobian!(u, U, F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    λ = x[end]
    evaluate_and_jacobian!(u, U, F.F, x, F.l, c.cache)
    u[end] = u[end] * λ + τ[1] * F.γ
    for j = 1:length(x)-1
        U[end, j] *= λ
    end
    U[end, end] = F.l[end]
    for i = 1:length(x)-1
        U[end, end] += F.l[i] * x[i]
    end
    nothing
end

function differentiate_parameters!(u, F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    u .= zero(eltype(u))
    u[end, 1] = F.γ
    u
end

function differentiate_parameters(F::TraceTestPencil, x, τ, c::TraceTestPencilCache)
    u = similar(c.cache.u, size(F)[1], 1)
    differentiate_parameters!(u, F, x, τ, c)
end
