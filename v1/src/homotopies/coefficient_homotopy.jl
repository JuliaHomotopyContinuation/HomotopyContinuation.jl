export CoefficientHomotopy, CoefficientHomotopyCache, gamma

@enum CoeffHomotopyActiveCoeffs begin
    CH_COEFFS_EVAL
    CH_COEFFS_DT
    CH_COEFFS_UNKNOWN
end


"""
    CoefficientHomotopy(support, start_coeffs, target_coeffs)

Construct the homotopy ``H(x, t) = ∑_{a ∈ Aᵢ} (c_a t + (1-t)d_a) x^a`` where ``c_a`` are
the start coefficients and ``d_a`` the target coefficients.
"""
struct CoefficientHomotopy{S<:AbstractSystem,T1,T2} <: AbstractHomotopy
    system::S
    start_coeffs::Vector{Vector{T1}}
    target_coeffs::Vector{Vector{T2}}
end

function CoefficientHomotopy(
    exponents_polys::Union{<:Vector{<:Matrix},<:MPPolys},
    start::Vector{<:Vector},
    target::Vector{<:Vector},
)
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

_nterms(A::Matrix) = size(A, 2)
_nterms(f::MPPoly) = MP.nterms(f)

(H::CoefficientHomotopy)(x, t, c = cache(H, x, t)) = evaluate(H, x, t, c)


"""
    CoefficientHomotopyCache

An simple cache for `StartTargetHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct CoefficientHomotopyCache{C<:AbstractSystemCache,T} <: AbstractHomotopyCache
    system::C
    coeffs::Vector{Vector{T}}
    diffs::Vector{Vector{T}}
    t::Base.RefValue{ComplexF64}
    active_coeffs::Base.RefValue{CoeffHomotopyActiveCoeffs}
end

function cache(H::CoefficientHomotopy{S,T1,T2}, x, t) where {S,T1,T2}
    system = cache(H.system, x, t)
    U = promote_type(T1, T2, typeof(t))
    coeffs = map(c -> convert.(U, c), H.start_coeffs)
    diffs = map(H.start_coeffs, H.target_coeffs) do start, target
        convert(Vector{U}, start - target)
    end
    active_coeffs = Ref(CH_COEFFS_UNKNOWN)
    CoefficientHomotopyCache(system, coeffs, diffs, Ref(complex(NaN)), active_coeffs)
end

Base.size(H::CoefficientHomotopy) = size(H.system)

function update_coeffs!(cache::CoefficientHomotopyCache, H::CoefficientHomotopy, t)
    if t == cache.t[]
        if cache.active_coeffs[] != CH_COEFFS_EVAL
            set_coefficients!(H.system, cache.coeffs)
            cache.active_coeffs[] = CH_COEFFS_EVAL
        end
        return nothing
    end

    if isreal(t)
        rt = real(t)
        rt1 = 1.0 - t
        @inbounds for k = 1:length(cache.coeffs)
            c = cache.coeffs[k]
            start = H.start_coeffs[k]
            target = H.target_coeffs[k]
            for i in eachindex(c)
                c[i] = rt * start[i] + rt1 * target[i]
            end
        end
    else
        t1 = 1.0 - t
        @inbounds for k = 1:length(cache.coeffs)
            c = cache.coeffs[k]
            start = H.start_coeffs[k]
            target = H.target_coeffs[k]
            for i in eachindex(c)
                c[i] = t * start[i] + t1 * target[i]
            end
        end
    end
    set_coefficients!(H.system, cache.coeffs)
    cache.active_coeffs[] = CH_COEFFS_EVAL
    cache.t[] = t
    nothing
end

function update_coeffs_dt!(cache::CoefficientHomotopyCache, H::CoefficientHomotopy, t)
    if cache.active_coeffs[] != CH_COEFFS_DT
        set_coefficients!(H.system, cache.diffs)
        cache.active_coeffs[] = CH_COEFFS_DT
    end
    nothing
end

function evaluate!(u, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function dt!(u, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs_dt!(c, H, t)
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function jacobian!(U, H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    @inbounds jacobian!(U, H.system, x, c.system)
end

function evaluate_and_jacobian!(
    u,
    U,
    H::CoefficientHomotopy,
    x,
    t,
    c::CoefficientHomotopyCache,
)
    update_coeffs!(c, H, t)
    @inbounds evaluate_and_jacobian!(u, U, H.system, x, c.system)
    nothing
end

function evaluate(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    evaluate(H.system, x, c.system)
end

function dt(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs_dt!(c, H, t)
    evaluate(H.system, x, c.system)
end

function jacobian(H::CoefficientHomotopy, x, t, c::CoefficientHomotopyCache)
    update_coeffs!(c, H, t)
    jacobian(H.system, x, c.system)
end
