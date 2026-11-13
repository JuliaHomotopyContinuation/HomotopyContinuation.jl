# The `scaled` harness mode: G(y) = c .* F(D y) with row factors cᵢ and column factors
# D = diag(dⱼ) from 10^U(-3, 3) with a fixed seed (see run.jl).
using Random

function harness_scaling(sys)
    rng = MersenneTwister(0x5ca1e)
    c = exp10.(6 .* rand(rng, length(sys.eqs)) .- 3)
    d = exp10.(6 .* rand(rng, length(sys.vars)) .- 3)
    c, d
end

# G(y; p) = c .* F(d .* y; p) around an unchanged system F (scaling only at the boundary).
# ComplexF64 inputs go through preallocated buffers so the scaling does not allocate on
# the tracker's hot path; other element types (the extended precision residual) allocate.
mutable struct Scaled{S<:AbstractSystem} <: AbstractSystem
    F::S
    c::Vector{Float64}
    d::Vector{Float64}
    x::Vector{ComplexF64}
    tx::Matrix{ComplexF64} # Taylor coefficients of x = d .* y, resized on demand
end
Scaled(F, c, d) =
    Scaled(F, c, d, zeros(ComplexF64, length(d)), zeros(ComplexF64, 0, length(d)))
Base.size(G::Scaled) = size(G.F)
ModelKit.variables(G::Scaled) = variables(G.F)
ModelKit.parameters(G::Scaled) = parameters(G.F)
# symbolic G, which `solve` builds to check for homogeneity
function ModelKit.System(G::Scaled)
    f = System(G.F)
    x, p = variables(f), parameters(f)
    System(G.c .* f(G.d .* x, p); variables = x, parameters = p)
end

unscale(G::Scaled, y::AbstractVector{ComplexF64}) = (G.x .= G.d .* y; G.x)
unscale(G::Scaled, y::AbstractVector) = G.d .* y
function unscale(G::Scaled, ty::TaylorVector{N,ComplexF64}) where {N}
    size(G.tx) == size(ty.data) || (G.tx = similar(ty.data))
    G.tx .= ty.data .* transpose(G.d)
    TaylorVector{N}(G.tx)
end
unscale(G::Scaled, ty::TaylorVector{N}) where {N} = TaylorVector{N}(ty.data .* transpose(G.d))

function ModelKit.evaluate!(u, G::Scaled, y, p = nothing)
    evaluate!(u, G.F, unscale(G, y), p)
    u .*= G.c
    u
end
function ModelKit.evaluate_and_jacobian!(u, U, G::Scaled, y, p = nothing)
    evaluate_and_jacobian!(u, U, G.F, unscale(G, y), p)
    u .*= G.c
    U .*= G.c .* transpose(G.d)
    nothing
end
function ModelKit.taylor!(u, v::Val, G::Scaled, ty, p = nothing)
    taylor!(u, v, G.F, unscale(G, ty), p)
    u .*= G.c
    u
end
