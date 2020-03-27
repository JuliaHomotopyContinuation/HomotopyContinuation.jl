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

function g!(u, H, x, t, (x¹,)::NTuple{1}, h, xh)
    xh .= x .+ h .* x¹
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, x, t, (ẋ, ẍ)::NTuple{2}, h, xh)
    h² = h^2
    @inbounds for i in eachindex(x)
        xh[i] = x[i] + h * ẋ[i] + h² * ẍ[i]
    end
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, x, t, (ẋ, ẍ, x³)::NTuple{3}, h, xh)
    h² = h^2
    h³ = h * h²
    @inbounds for i in eachindex(x)
        xh[i] = x[i] + h * ẋ[i] + h² * ẍ[i] + h³ * x³[i]
    end
    evaluate!(u, H, xh, t + h)
end

function diff_t!(
    u,
    H::AbstractHomotopy,
    x,
    t,
    dx::Tuple{},
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_3 = 6.055454452393343e-6
    h_1_5 = 0.000740095979741405

    if 2τ > h_1_5
        # apply S^{1,1,4} formula
        h = h_1_5
        evaluate!(u₁, H, x, t + h)
        evaluate!(u₂, H, x, t - h)
        evaluate!(u₃, H, x, t + 2h)
        evaluate!(u₄, H, x, t - 2h)
        u .= (2.0 .* (u₁ .- u₂) ./ 3.0 .- (u₃ .- u₄) ./ 12.0) ./ h
    elseif τ > h_1_3  || !use_extended_precision
        # apply S^{1,1,2} formula
        h = min(τ, h_1_3)
        evaluate!(u₁, H, x, t + h)
        evaluate!(u₂, H, x, t - h)
        u .= 0.5 .* (u₁ .- u₂) ./ h
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{1,1,2} formula
        h = min(τ, 2.3099725541661633e-11) #eps(DoubleF64)^(1/3)
        evaluate!(u₁_extended, H, xh_extended, ComplexDF64(t + h))
        evaluate!(u₂_extended, H, xh_extended, ComplexDF64(t - h))
        u .= 0.5 .* (u₁_extended .- u₂_extended) ./ h
    end

    u
end

function diff_t!(
    u,
    H::AbstractHomotopy,
    x,
    t,
    dx::NTuple{1},
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_4 = 0.0001220703125 # eps()^(1/(2+2))
    h_1_6 = 0.002460783300575925 # eps()^(1/(2+4))

    if 0.5τ > h_1_6
        h = h_1_6
        # apply S^{2,2,4} formula
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        g!(u₃, H, x, t, dx, 2h, xh)
        g!(u₄, H, x, t, dx, -2h, xh)
        h2 = h^2
        u .= (2.0 .* (u₁ .+ u₂) ./ 3.0 .- (u₃ .+ u₄) ./ 24.0) ./ h2
    elseif τ > h_1_4  || !use_extended_precision
        h = min(τ, h_1_4)
        # apply S^{2,2,2} formula
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        h2 = h^2
        u .= 0.5 .* (u₁ .+ u₂) ./ h2
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{2,2,2} formula
        h = min(τ, 1.0536712127723509e-8) #eps(DoubleF64)^(1/4)
        g!(u₁_extended, H, x, t, dx, h, xh_extended)
        g!(u₂_extended, H, x, t, dx, -h, xh_extended)
        h2 = h^2
        u .= 0.5 .* (u₁_extended .+ u₂_extended) ./ h2
    end

    u
end

function diff_t!(
    u,
    H::AbstractHomotopy,
    x,
    t,
    dx::NTuple{2},
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_5 = 0.000740095979741405 # eps^(1/(3+2))
    h_1_7 = 0.0058046651919412065 # eps()^(1/(3+4))
    if 0.5τ > h_1_7
        # apply S^{3,3,4} formula
        h = h_1_7
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        g!(u₃, H, x, t, dx, 2h, xh)
        g!(u₄, H, x, t, dx, -2h, xh)
        # divide by 6 since we compute taylor series
        h3 = 6 * h^3
        u .= (4.0 .* (u₁ .- u₂) .- 0.125 .* (u₃ .- u₄)) ./ h3
    elseif τ > h_1_5 || !use_extended_precision
        h = min(τ, h_1_5)
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        h3 = h^3
        u .= 0.5 .* (u₁ .- u₂) ./ h3
    else
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        #use extended precision
        # apply S^{3,3,2} formula
        h = min(τ, 4.151108566742532e-7) #eps(DoubleF64)^(1/5)
        g!(u₁_extended, H, x, t, dx, h, xh_extended)
        g!(u₂_extended, H, x, t, dx, -h, xh_extended)
        h3 = h^3
        u .= 0.5 .* (u₁_extended .- u₂_extended) ./ h3
    end

    u
end

function diff_t!(
    u,
    H::AbstractHomotopy,
    x,
    t,
    dx::NTuple{3},
    ND::NumericalDifferentiation,
    τ::Float64 = Inf,
    use_extended_precision::Bool = true
)
    @unpack u₁, u₂, u₃, u₄, xh = ND

    h_1_6 = 0.002460783300575925
    h_1_8 = 0.011048543456039806
    if 0.5τ > h_1_8
        # apply S^{4,4,4} formula
        h = h_1_8
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        g!(u₃, H, x, t, dx, 2h, xh)
        g!(u₄, H, x, t, dx, -2h, xh)
        h2 = h^2
        h4 = 24 * h2^2
        # divide by 4! since we compute taylor series
        u .= (16.0 .* (u₁ .+ u₂) .- 0.25 .* (u₃ .+ u₄)) ./ h4
    elseif τ > h_1_6 || !use_extended_precision
        h = min(τ, h_1_6)
        g!(u₁, H, x, t, dx, h, xh)
        g!(u₂, H, x, t, dx, -h, xh)
        h2 = h^2
        h4 = h2^2
        u .= 0.5 .* (u₁ .+ u₂) ./ h4
    else
        #use extended precision
        @unpack u₁_extended, u₂_extended, xh_extended = ND
        # apply S^{4,4,2} formula
        h = min(τ, 4.806217383937355e-6) #eps(DoubleF64)^(1/6)
        g!(u₁_extended, H, x, t, dx, h, xh_extended)
        g!(u₂_extended, H, x, t, dx, -h, xh_extended)
        h2 = h^2
        h4 = h2^2
        u .= 0.5 .* (u₁_extended .+ u₂_extended) ./ h4
    end
    u
end

## Default handling
diff_t!(u, H::AbstractHomotopy, x, t, dx::Tuple, DS::AutomaticDifferentiation, τ) =
    diff_t!(u, H, x, t, dx)

@generated function diff_t!(u, H, x, t, dx::NTuple{M}, AD::Val{N}, ND, τ; use_extended_precision::Bool = true) where {M,N}
    if M < N
        :(diff_t!(u, H, x, t, dx, AutomaticDifferentiation(), τ))
    else
        :(diff_t!(u, H, x, t, dx, ND, τ, use_extended_precision))
    end
end
