export PolyhedralHomotopy, PolyhedralHomotopyCache, gamma

"""
    PolyhedralHomotopy(G, F; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, s) = Σ_α c_α e^{(α̂⋅γ̂) s} x^α``.
"""
mutable struct PolyhedralHomotopy{S<:SPSystem, T} <: AbstractHomotopy
    system::S
    support::Vector{Matrix{Int32}}
    lifting::Vector{Vector{Int32}}
    coeffs::Vector{Vector{T}}
    s_weights::Vector{Vector{Float64}}
    π::Vector{Vector{Int}}
end

function PolyhedralHomotopy(support::Vector{<:Matrix}, lifting::Vector{Vector{Int32}}, coeffs::Vector{Vector{T}}) where T
    @assert all(length.(coeffs) .== size.(support,2))
    @assert length(lifting) == length(support)

    system = SPSystem(support, coeffs)
    # StaticPolynomials possibly changes the order of the terms.
    # Therefore the coefficient vectors are maybe no more correct.
    # We correct this by applying the permutation applied to the columns
    # of the support
    coeffs₀ = deepcopy(coeffs)
    s_weights = [zeros(Float64, length(c)) for c in coeffs]
    π = Vector{Vector{Int}}()
    for (i, f) in enumerate(system.system.polys)
        p = SP.permutation(f)
        push!(π, p)
        permute!(coeffs₀[i], p)
    end

    PolyhedralHomotopy(system, convert(Vector{Matrix{Int32}}, support), lifting, coeffs₀, s_weights, π)
end

function update_cell!(H::PolyhedralHomotopy, cell::MixedSubdivisions.MixedCell)
    n = length(H.lifting)
    γ = cell.normal
    for i in 1:n
        ω = H.lifting[i]
        βᵢ = cell.β[i]
        Aᵢ = H.support[i]
        m = size(Aᵢ, 2)
        for j = 1:m
            sᵢⱼ = Float64(ω[j])
            for k = 1:n
                sᵢⱼ += Aᵢ[k,j] * γ[k]
            end
            sᵢⱼ -= βᵢ
            H.s_weights[i][j] = sᵢⱼ
        end
        permute!(H.s_weights[i], H.π[i])
    end
    H
end

"""
    PolyhedralHomotopyCache

An simple cache for `PolyhedralHomotopyCache`.
"""
mutable struct PolyhedralHomotopyCache{C<:AbstractSystemCache, T} <: AbstractHomotopyCache
    system::C
    coeffs::Vector{Vector{T}}
    coeffs_dt::Vector{Vector{T}}
    s::Float64
    ds::Float64
    active_coeffs::ActiveCoeffs
end

function cache(H::PolyhedralHomotopy{S,T}, x, s) where {S,T}
    system = cache(H.system, x, s)
    U = promote_type(T, typeof(s))
    coeffs = map(c -> convert.(U, c), H.coeffs)
    s = ds = NaN
    active_coeffs = COEFFS_UNKNOWN
    PolyhedralHomotopyCache(system, coeffs, deepcopy(coeffs), s, ds, active_coeffs)
end

Base.size(H::PolyhedralHomotopy) = size(H.system)

function update_coeffs!(cache::PolyhedralHomotopyCache, H::PolyhedralHomotopy, s)
    if s == cache.s
        if cache.active_coeffs != COEFFS_EVAL
            set_coefficients!(H.system, cache.coeffs)
            cache.active_coeffs == COEFFS_EVAL
        end
        return nothing
    end

    cache.s = s
    @inbounds for k in 1:length(cache.coeffs)
        e = cache.coeffs[k]
        for i in eachindex(e)
            e[i] = exp(H.s_weights[k][i] * s) * H.coeffs[k][i]
        end
    end
    set_coefficients!(H.system, cache.coeffs)
    cache.active_coeffs == COEFFS_EVAL
    nothing
end

function update_coeffs_dt!(cache::PolyhedralHomotopyCache, H::PolyhedralHomotopy, s)
    if s == cache.ds
        if cache.active_coeffs != COEFFS_DT
            set_coefficients!(H.system, cache.coeffs)
            cache.active_coeffs == COEFFS_DT
        end
        return nothing
    end

    if s == cache.s
        cache.ds = s
        @inbounds for k in 1:length(cache.coeffs_dt)
            e = cache.coeffs_dt[k]
            for i in eachindex(e)
                e[i] = H.s_weights[k][i] * cache.coeffs[k][i]
            end
        end
    else
        cache.ds = s
        @inbounds for k in 1:length(cache.coeffs_dt)
            e = cache.coeffs_dt[k]
            for i in eachindex(e)
                e[i] = H.s_weights[k][i] * exp(H.s_weights[k][i] * s) * H.coeffs[k][i]
            end
        end
    end
    set_coefficients!(H.system, cache.coeffs_dt)
    cache.active_coeffs == COEFFS_DT
    nothing
end


function evaluate!(u, H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs!(c, H, real(s))
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function dt!(u, H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs_dt!(c, H, real(s))
    @inbounds evaluate!(u, H.system, x, c.system)
    u
end

function jacobian!(U, H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs!(c, H, real(s))
    @inbounds jacobian!(U, H.system, x, c.system)
end

function evaluate_and_jacobian!(u, U, H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs!(c, H, real(s))
    @inbounds evaluate_and_jacobian!(u, U, H.system, x, c.system)
    nothing
end

function evaluate(H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs!(c, H, real(s))
    evaluate(H.system, x, c.system)
end

function dt(H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs_dt!(c, H, real(s))
    evaluate(H.system, x, c.system)
end

function jacobian(H::PolyhedralHomotopy, x, s, c::PolyhedralHomotopyCache)
    update_coeffs!(c, H, real(s))
    jacobian(H.system, x, c.system)
end

(H::PolyhedralHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)
