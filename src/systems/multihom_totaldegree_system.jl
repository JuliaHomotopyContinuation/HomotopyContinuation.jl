export MultiHomTotalDegreeSystem

"""
    MultiHomogeneousTotalDegreeSystem(polynomials, vars) <: AbstractSystem

Create a multi-homogeneous total degree system as described in

An efficient start system for multi-homogeneous polynomial continuation,
Wampler, C.W. Numer. Math. (1993) 66: 517. https://doi.org/10.1007/BF01385710
"""
struct MultiHomTotalDegreeSystem{T} <: AbstractSystem
    D::Matrix{Int}
    C::Matrix{Float64}
    scaling_factors::Vector{T}
end

function MultiHomTotalDegreeSystem(D, C)
    MultiHomTotalDegreeSystem(D, C, ones(size(D, 2)))
end

struct MultiHomTotalDegreeSystemCache{M,T} <: AbstractSystemCache
    B::Matrix{T}
    Ĝ::Matrix{T}
    R::Matrix{T}
    S::Matrix{T}
    ranges::NTuple{M,UnitRange{Int}}
    homvars::NTuple{M,Int}
end

Base.size(F::MultiHomTotalDegreeSystem) = (size(F.D, 2), sum(size(F.D)))
Base.length(F::MultiHomTotalDegreeSystem) = size(F.D, 2)

function cache(F::MultiHomTotalDegreeSystem, x::PVector{TX,M}) where {TX,M}
    T = typeof(F.C[1, 1] * x[1] + F.C[1, 1] * x[1])
    n = length(F)
    B, Ĝ = zeros(T, M, n), zeros(T, M, n)
    R, S = zeros(T, M, n), zeros(T, M, n)
    ranges_homvars = ProjectiveVectors.dimension_indices_homvars(x)
    ranges = first.(ranges_homvars)
    homvars = last.(ranges_homvars)
    MultiHomTotalDegreeSystemCache(B, Ĝ, R, S, ranges, homvars)
end

function evaluate!(
    u,
    F::MultiHomTotalDegreeSystem,
    z::ProjectiveVectors.PVector,
    cache::MultiHomTotalDegreeSystemCache{M},
) where {M}
    D, C = F.D, F.C
    B, ranges, homvars = cache.B, cache.ranges, cache.homvars

    n = size(B, 2)

    @boundscheck checkbounds(u, 1:n)

    # Compute all bᵢⱼ and store in B
    # Since B is column major we store bᵢⱼ in B[j, i]
    @inbounds for i = 1:n
        for j = 1:M
            if D[j, i] != 0
                bᵢⱼ = -zero(B[j, i])
                for k in ranges[j]
                    bᵢⱼ += C[k, i] * z[k]
                end
                B[j, i] = bᵢⱼ
            end
        end
    end

    @inbounds for i = 1:n
        gᵢ = one(eltype(u))
        for j = 1:M
            dᵢⱼ, bᵢⱼ = D[j, i], B[j, i]
            if dᵢⱼ != 0
                gᵢ *= bᵢⱼ^dᵢⱼ - z[homvars[j]]^dᵢⱼ
            end
        end
        u[i] = gᵢ * F.scaling_factors[i]
    end

    u
end
function evaluate(F::MultiHomTotalDegreeSystem, x, cache::MultiHomTotalDegreeSystemCache)
    evaluate!(similar(x, size(F, 1)), F, x, cache)
end

function jacobian!(
    U,
    F::MultiHomTotalDegreeSystem,
    z,
    cache::MultiHomTotalDegreeSystemCache{M},
) where {M}
    evaluate_and_jacobian!(nothing, U, F, z, cache)
    U
end
function jacobian(F::MultiHomTotalDegreeSystem, x, cache::MultiHomTotalDegreeSystemCache)
    jacobian!(similar(x, size(F)), F, x, cache)
end

function evaluate_and_jacobian!(
    u,
    U,
    F::MultiHomTotalDegreeSystem,
    z,
    cache::MultiHomTotalDegreeSystemCache{M},
) where {M}
    n, N = size(F)

    @boundscheck checkbounds(U, 1:n, 1:N)
    if u !== nothing
        @boundscheck checkbounds(u, 1:n)
    end

    D, C = F.D, F.C
    B, Ĝ, R, S, ranges, homvars = cache.B,
        cache.Ĝ,
        cache.R,
        cache.S,
        cache.ranges,
        cache.homvars

    n = size(B, 2)

    @inbounds for i = 1:n
        # Compute all bᵢⱼ and store in B
        # Since B is column major we store bᵢⱼ in B[j, i]
        for j = 1:M
            if D[j, i] != 0
                bᵢⱼ = -zero(B[j, i])
                for k in ranges[j]
                    bᵢⱼ += C[k, i] * z[k]
                end
                B[j, i] = bᵢⱼ
            end
        end

        # Compute all ĝᵢⱼ and store in Ĝ
        # Since Ĝ is column major we store ĝᵢⱼ in Ĝ[j, i]
        for j = 1:M
            dᵢⱼ, bᵢⱼ = D[j, i], B[j, i]
            if dᵢⱼ == 0
                Ĝ[j, i] = one(eltype(Ĝ))
            else
                Ĝ[j, i] = bᵢⱼ^dᵢⱼ - z[homvars[j]]^dᵢⱼ
            end
        end

        # Accumulate subproducts forward
        R[1, i] = rᵢⱼ_prev = Ĝ[1, i]
        for j = 2:M
            if D[j, i] != 0 # otherwise Ĝ[j, i] = 1
                R[j, i] = rᵢⱼ_prev = rᵢⱼ_prev * Ĝ[j, i]
            else
                R[j, i] = rᵢⱼ_prev
            end
            if u !== nothing
                u[i] = rᵢⱼ_prev * F.scaling_factors[i]
            end
        end

        # Accumulate subproducts backward
        S[M, i] = sᵢⱼ_prev = Ĝ[M, i]
        for j = M-1:-1:1
            if D[j, i] != 0 # otherwise Ĝ[j, i] = 1
                S[j, i] = sᵢⱼ_prev = sᵢⱼ_prev * Ĝ[j, i]
            else
                S[j, i] = sᵢⱼ_prev
            end
        end

        # Compute partial derivatives
        for j = 1:M
            dᵢⱼ = D[j, i]
            for k in ranges[j]
                c = C[k, i]
                if iszero(c) || iszero(dᵢⱼ)
                    U[i, k] = zero(eltype(U))
                else
                    if dᵢⱼ == 1
                        u_ik = c
                    else
                        u_ik = dᵢⱼ * B[j, i]^(dᵢⱼ - 1) * c
                    end
                    if j > 1
                        u_ik *= R[j-1, i]
                    end
                    if j < M
                        u_ik *= S[j+1, i]
                    end
                    U[i, k] = u_ik * F.scaling_factors[i]
                end
            end

            k = homvars[j]
            if iszero(dᵢⱼ)
                U[i, k] = zero(eltype(U))
            else
                if dᵢⱼ == 1
                    u_ik = -dᵢⱼ
                else
                    u_ik = -dᵢⱼ * z[k]^(dᵢⱼ - 1)
                end
                if j > 1
                    u_ik *= R[j-1, i]
                end
                if j < M
                    u_ik *= S[j+1, i]
                end
                U[i, k] = u_ik * F.scaling_factors[i]
            end
        end
    end
    nothing
end
