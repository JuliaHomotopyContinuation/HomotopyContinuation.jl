export CoefficientHomotopy, CoefficientHomotopyCache, gamma

"""
    CoefficientHomotopy(G, F; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct CoefficientHomotopy{S<:AbstractSystem, T} <: AbstractHomotopy
    system::S
    start_coeffs::Vector{Vector{T}}
    target_coeffs::Vector{Vector{T}}
end

function CoefficientHomotopy(exponents_polys::Union{<:Vector{<:Matrix}, <:MPPolys}, start::Vector{Vector{T}}, target::Vector{Vector{T}}) where T
    @assert all(length.(start) .== length.(target) .== _nterms.(exponents_polys))

    system = SPSystem(exponents_polys, start)
    # StaticPolynomials possibly changes the order of the terms.
    # Therefore the coefficient vectors are maybe no more correct.
    # We correct this by applying the permutation applied to the columns
    # of the support
    start_coeffs = deepcopy(start)
    target_coeffs = deepcopy(target)
    for (i, f) in enumerate(system.system.polys)
        p = SP.permutation(f)
        permute!(start_coeffs[i], p)
        permute!(target_coeffs[i], p)
    end

    CoefficientHomotopy(system, start_coeffs, target_coeffs)
end

_nterms(A::Matrix) = size(A,2)
_nterms(f::MPPoly) = MP.nterms(f)

(H::CoefficientHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

"""
    CoefficientHomotopyCache

An simple cache for `StartTargetHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct CoefficientHomotopyCache{C<:AbstractSystemCache, T} <: AbstractHomotopyCache
    system::C
    coeffs::Vector{Vector{T}}
    diffs::Vector{Vector{T}}
    t::ComplexF64
end

function cache(H::CoefficientHomotopy{S,T}, x, t) where {S,T}
    system = cache(H.system, x, t)
    U = promote_type(T, typeof(t))
    coeffs = map(c -> convert.(U, c), H.start_coeffs)
    diffs = map(H.start_coeffs, H.target_coeffs) do start, target
        convert(Vector{U}, start - target)
    end
    CoefficientHomotopyCache(system, coeffs, diffs, complex(1.0))
end

Base.size(H::CoefficientHomotopy) = size(H.system)

function update_coeffs!(cache::CoefficientHomotopyCache, H::CoefficientHomotopy, t)
    t == cache.t && return nothing

    for k in 1:length(cache.coeffs)
        c = cache.coeffs[k]
        start = H.start_coeffs[k]
        target = H.target_coeffs[k]
        for i in eachindex(c)
            c[i] = t * start[i] + (1.0 - t) * target[i]
        end
    end
    set_coefficients!(H.system, cache.coeffs)
    nothing
end

function evaluate!(u, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function dt!(u, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    set_coefficients!(H.system, c.diffs)
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function jacobian!(U, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    @inbounds jacobian!(U, H.system, x, c.system)
end

function evaluate_and_jacobian!(u, U, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    @inbounds evaluate_and_jacobian!(u, U, H.system, x, c.system)
    nothing
end

function evaluate(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    evaluate(H.system, x, c.system)
end

function dt(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    set_coefficients!(H.system, c.diffs)
    evaluate(H.system, x, c.system)
end

function jacobian(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    jacobian(H.system, x, c.system)
end
