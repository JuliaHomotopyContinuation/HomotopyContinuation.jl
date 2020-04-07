abstract type DifferentiationStrategy end

struct AutomaticDifferentiation <: DifferentiationStrategy end

# Numerical Differentiation

struct NumericalDifferentiation{V<:AbstractVector,V̄<:AbstractVector} <:
       DifferentiationStrategy
    xh::V
    u₁::Vector{ComplexF64}
    u₂::Vector{ComplexF64}
    u₃::Vector{ComplexF64}
    u₄::Vector{ComplexF64}
    xh_extended::V̄
    u₁_extended::Vector{ComplexDF64}
    u₂_extended::Vector{ComplexDF64}
end
NumericalDifferentiation(x::AbstractVector, n::Int) = NumericalDifferentiation(
    copy(x),
    (zeros(ComplexF64, n) for i = 1:4)...,
    ComplexDF64.(x),
    (zeros(ComplexDF64, n) for i = 1:2)...,
)

function g!(u, H, tx::TaylorVector{1}, t, h, xh)
    for i = 1:length(tx)
        x, = tx[i]
        xh[i] = x
    end
    if t isa Complex
        evaluate!(u, H, xh, ComplexF64(t + h))
    else
        evaluate!(u, H, xh, Float64(t + h))
    end

end
function g!(u, H, tx::TaylorVector{2}, t, h, xh)
    for i = 1:length(tx)
        x, x¹ = tx[i]
        xh[i] = x
        xh[i] += h * x¹
    end
    if t isa Complex
        evaluate!(u, H, xh, ComplexF64(t + h))
    else
        evaluate!(u, H, xh, Float64(t + h))
    end
end
function g!(u, H, tx::TaylorVector{3}, t, h, xh)
    h² = h^2
    for i = 1:length(tx)
        x, x¹, x² = tx[i]
        xh[i] = x
        xh[i] += h * x¹
        xh[i] += h² * x²
    end
    if t isa Complex
        evaluate!(u, H, xh, ComplexF64(t + h))
    else
        evaluate!(u, H, xh, Float64(t + h))
    end
end

function g!(u, H, tx::TaylorVector{4}, t, h, xh)
    h² = h^2
    h³ = h * h²
    for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        xh[i] = x
        xh[i] += h * x¹
        xh[i] += h² * x²
        xh[i] += h³ * x³
    end
    if t isa Complex
        evaluate!(u, H, xh, ComplexF64(t + h))
    else
        evaluate!(u, H, xh, Float64(t + h))
    end
end

function taylor!(
    u,
    ::Val{1},
    H::AbstractHomotopy,
    x::TaylorVector{1},
    t,
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true,
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_3 = 6.055454452393343e-6
    h_1_5 = 0.000740095979741405

    # @show τ
    if 2τ > h_1_5
        # apply S^{1,1,4} formula
        h = h_1_5
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        g!(u₃, H, x, t, 2h, xh)
        g!(u₄, H, x, t, -2h, xh)
        u .= (2.0 .* (u₁ .- u₂) ./ 3.0 .- (u₃ .- u₄) ./ 12.0) ./ h
    elseif τ > h_1_3 || !use_extended_precision
        # apply S^{1,1,2} formula
        h = min(τ, h_1_3)
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        u .= 0.5 .* (u₁ .- u₂) ./ h
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{1,1,2} formula
        h = DoubleF64(min(τ, 2.3099725541661633e-11)) #eps(DoubleF64)^(1/3)
        g!(u₁_extended, H, x, t, h, xh_extended)
        g!(u₂_extended, H, x, t, h, xh_extended)
        u .= 0.5 .* (u₁_extended .- u₂_extended) ./ h
    end

    u
end

function taylor!(
    u,
    ::Val{2},
    H::AbstractHomotopy,
    x::TaylorVector{2},
    t,
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true,
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_4 = 0.0001220703125 # eps()^(1/(2+2))
    h_1_6 = 0.002460783300575925 # eps()^(1/(2+4))

    if 0.5τ > h_1_6
        h = h_1_6
        # apply S^{2,2,4} formula
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        g!(u₃, H, x, t, 2h, xh)
        g!(u₄, H, x, t, -2h, xh)
        h2 = h^2
        u .= (2.0 .* (u₁ .+ u₂) ./ 3.0 .- (u₃ .+ u₄) ./ 24.0) ./ h2
    elseif τ > h_1_4 || !use_extended_precision
        h = min(τ, h_1_4)
        # apply S^{2,2,2} formula
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        h2 = h^2
        u .= 0.5 .* (u₁ .+ u₂) ./ h2
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{2,2,2} formula
        h = DoubleF64(min(τ, 1.0536712127723509e-8)) #eps(DoubleF64)^(1/4)
        g!(u₁_extended, H, x, t, h, xh_extended)
        g!(u₂_extended, H, x, t, -h, xh_extended)
        h2 = h^2
        u .= 0.5 .* (u₁_extended .+ u₂_extended) ./ h2
    end

    u
end

function taylor!(
    u,
    ::Val{3},
    H::AbstractHomotopy,
    x::TaylorVector{3},
    t,
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true,
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_5 = 0.000740095979741405 # eps^(1/(3+2))
    h_1_7 = 0.0058046651919412065 # eps()^(1/(3+4))
    if 0.5τ > h_1_7
        # apply S^{3,3,4} formula
        h = h_1_7
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        g!(u₃, H, x, t, 2h, xh)
        g!(u₄, H, x, t, -2h, xh)
        # divide by 6 since we compute taylor series
        h3 = 6 * h^3
        u .= (4.0 .* (u₁ .- u₂) .- 0.125 .* (u₃ .- u₄)) ./ h3
    elseif τ > h_1_5 || !use_extended_precision
        h = min(τ, h_1_5)
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        h3 = h^3
        u .= 0.5 .* (u₁ .- u₂) ./ h3
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{3,3,2} formula
        h = DoubleF64(min(τ, 4.151108566742532e-7)) #eps(DoubleF64)^(1/5)
        g!(u₁_extended, H, x, t, h, xh_extended)
        g!(u₂_extended, H, x, t, -h, xh_extended)
        h3 = h^3
        u .= 0.5 .* (u₁_extended .- u₂_extended) ./ h3
    end

    u
end

function taylor!(
    u,
    ::Val{4},
    H::AbstractHomotopy,
    x::TaylorVector{4},
    t,
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = false,
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_6 = 0.002460783300575925
    if τ > h_1_6 || !use_extended_precision
        h = min(τ, h_1_6)
        g!(u₁, H, x, t, h, xh)
        g!(u₂, H, x, t, -h, xh)
        h2 = h^2
        h4 = h2^2
        u .= 0.5 .* (u₁ .+ u₂) ./ h4
    else
        #use extended precision
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        # apply S^{4,4,2} formula
        h = DoubleF64(min(τ, 4.806217383937355e-6)) #eps(DoubleF64)^(1/6)
        g!(u₁_extended, H, x, t, h, xh_extended)
        g!(u₂_extended, H, x, t, -h, xh_extended)
        h2 = h^2
        h4 = h2^2
        u .= 0.5 .* (u₁_extended .+ u₂_extended) ./ h4
    end
    u
end

## Default handling ignores incremental
taylor!(u, v::Val, H::AbstractHomotopy, tx::TaylorVector, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)

@generated function taylor!(
    u,
    v::Val{M},
    H,
    tx,
    t,
    AD::Val{N},
    ND,
    τ;
    use_extended_precision::Bool = true,
    incremental::Bool = false,
) where {M,N}
    if M ≤ N
        :(taylor!(u, v, H, tx, t, incremental))
    else
        :(taylor!(u, v, H, tx, t, ND, τ, use_extended_precision))
    end
end
