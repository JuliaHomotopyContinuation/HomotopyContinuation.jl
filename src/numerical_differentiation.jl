mutable struct NumericalDifferentiation
    xh::Vector{ComplexF64}
    u₁::Vector{ComplexF64}
    u₂::Vector{ComplexF64}
    u₃::Vector{ComplexF64}
    h::Vector{Float64}
    default_h::Vector{Float64}
    Δ::NTuple{4,Vector{Float64}}
    tmp_Δ::Vector{Float64}
end

function NumericalDifferentiation(x::AbstractVector, n::Int)
    default_h = [eps()^(1 / (2 + k)) for k = 1:4]
    h = copy(default_h)
    NumericalDifferentiation(
        copy(x),
        (zeros(ComplexF64, n) for i = 1:3)...,
        h,
        default_h,
        ntuple(_ -> zeros(n), 4),
        zeros(n),
    )
end

function g!(u, H, tx::TaylorVector{1}, t, h, xh)
    @inbounds for i = 1:length(tx)
        xh[i] = first(tx[i])
    end
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, tx::TaylorVector{2}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹ = tx[i]
        xh[i] = muladd(h, x¹, x)
    end
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, tx::TaylorVector{3}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹, x² = tx[i]
        xh[i] = muladd(h, muladd(h, x², x¹), x)
    end
    evaluate!(u, H, xh, t + h)
end

function g!(u, H, tx::TaylorVector{4}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        xh[i] = muladd(h, muladd(h, muladd(h, x³, x²), x¹), x)
    end
    evaluate!(u, H, xh, t + h)
end

function finite_diff!(
    u,
    Δ,
    H::AbstractHomotopy,
    x::TaylorVector,
    t,
    ND::NumericalDifferentiation;
    order::Int,
    dist_to_target::Float64,
)
    h̄ = ND.h[order]
    if dist_to_target < h̄
        finite_diff!(u, Δ, H, x, t, ND, im * h̄; order = order)
    else
        finite_diff!(u, Δ, H, x, t, ND, h̄; order = order)
    end
end
function finite_diff!(
    u,
    Δ,
    H::AbstractHomotopy,
    x::TaylorVector,
    t,
    ND::NumericalDifferentiation,
    h;
    order::Int,
)
    @unpack u₁, u₂, xh = ND
    N = order

    g!(u₁, H, x, t, h, ND.xh)
    g!(u₂, H, x, t, -h, ND.xh)

    hN = h^N
    if iseven(N)
        u .= 0.5 .* (u₁ .+ u₂) ./ hN
    else
        u .= 0.5 .* (u₁ .- u₂) ./ hN
    end

    δ = 0.0
    ntrusted = 0
    abs_hN = fast_abs(hN)
    @inbounds for i in eachindex(u)
        fwdᵢ = fast_abs(u₁[i]) / abs_hN
        bwdᵢ = fast_abs(u₂[i]) / abs_hN
        cᵢ = fast_abs(u[i])
        fδᵢ = abs(1.0 - fwdᵢ / cᵢ)
        bδᵢ = abs(1.0 - bwdᵢ / cᵢ)

        if !iszero(fwdᵢ) && !iszero(bwdᵢ) && !iszero(cᵢ)
            Δ[i] = max(fδᵢ, bδᵢ)
        else
            Δ[i] = 0.0
        end
    end
    δ
end

function best_h_finite_diff!(
    u,
    Δ,
    H::AbstractHomotopy,
    x::TaylorVector,
    t,
    ND::NumericalDifferentiation;
    order::Int,
    dist_to_target::Float64,
)
    h = ND.default_h[order]
    max_δ = -Inf
    min_h = h
    max_ntrusted = 0
    for k = 0:2
        ND.h[order] = h
        finite_diff!(
            ND.u₃,
            ND.tmp_Δ,
            H,
            x,
            t,
            ND;
            order = order,
            dist_to_target = dist_to_target,
        )
        δ̂ = sum(exp ∘ -, ND.tmp_Δ) / size(H, 1)
        if δ̂ > max_δ
            max_δ = δ̂
            min_h = h
            u .= ND.u₃
            Δ .= ND.tmp_Δ
        end
        h *= 0.0625
    end
    ND.h[order] = min_h
    h = ND.h[order]
end

function taylor!(
    u,
    ::Val{N},
    H::AbstractHomotopy,
    x::TaylorVector{N},
    t,
    ND::NumericalDifferentiation;
    dist_to_target::Float64,
) where {N}
    finite_diff!(u, ND.Δ[N], H, x, t, ND; order = N, dist_to_target = dist_to_target)
    tol = min(0.1, 10 * ND.default_h[N])
    δ = maximum(ND.Δ[N])
    if δ > tol
        h = ND.h[N]
        best_h_finite_diff!(
            u,
            ND.Δ[N],
            H,
            x,
            t,
            ND;
            order = N,
            dist_to_target = dist_to_target,
        )
        h = ND.h[N]
        return maximum(ND.Δ[N])
    else
        return δ
    end
end

## Default handling ignores incremental
taylor!(u, v::Val, H::AbstractHomotopy, tx::TaylorVector, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)

## Type dispatch on automatic differentiation or numerical differentiation
@generated function taylor!(
    u,
    v::Val{M},
    H,
    tx,
    t,
    AD::Val{N},
    ND::NumericalDifferentiation;
    cond::Float64,
    dist_to_target::Float64,
    incremental::Bool = false,
) where {M,N}
    if M ≤ N
        quote
            taylor!(u, v, H, tx, t, incremental)
            true
        end
    else
        quote
            δ = taylor!(u, v, H, tx, t, ND; dist_to_target = dist_to_target)
            δ * cond < 1
        end
    end
end
