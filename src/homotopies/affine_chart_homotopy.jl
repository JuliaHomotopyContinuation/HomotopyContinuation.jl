export AffineChartHomotopy, on_affine_chart

"""
    AffineChartHomotopy(H::AbstractHomotopy, v::PVector{T,N})

Given a homotopy ``H(x,t): (ℙ^{m_1} × ⋯ × ℙ^{m_N}) × ℂ → ℂⁿ`` this creates a new affine
homotopy ``H̄`` which operates on the affine chart defined by the vector
``v ∈ ℙ^{m_1} × ⋯ × ℙ^{m_N}`` and the augmented conditions ``vᵀx = 1``.
"""
struct AffineChartHomotopy{H<:AbstractHomotopy,N} <: AbstractHomotopy
    homotopy::H
    chart::PVector{ComplexF64,N}
end

"""
    on_affine_chart(H::AbstractHomotopy, proj_dims)

Construct an `AffineChartHomotopy` on a randomly generated chart `v`. Each entry is drawn
idepdently from a univariate normal distribution.
"""
function on_affine_chart(H::AbstractHomotopy, proj_dims = [size(H,2) - 1])
    N = length(proj_dims)
    chart = PVector(randn(ComplexF64, sum(proj_dims) + N), tuple(proj_dims...))
    AffineChartHomotopy(H, chart)
end

function Base.size(H::AffineChartHomotopy{<:AbstractHomotopy,N}) where {N}
    m, n = size(H.homotopy)
    (m + N, n)
end

function set_solution!(x::AbstractVector, H::AffineChartHomotopy, y::AbstractVector, t)
    x .= y
    on_chart!(x, H.chart)
    x
end
function on_chart!(x::AbstractVector, v::PVector)
    ranges = dimension_indices(v)
    for range in ranges
        λ = zero(eltype(x))
        @inbounds for i in range
            λ += v[i] * x[i]
        end
        λ⁻¹ = @fastmath inv(λ)
        for i in range
            x[i] *= λ⁻¹
        end
    end
    x
end

function evaluate!(u, H::AffineChartHomotopy{<:Any,N}, x::AbstractVector, t) where {N}
    evaluate!(u, H.homotopy, x, t)
    evaluate_chart!(u, H.chart, x)
    u
end

function evaluate_and_jacobian!(
    u,
    U,
    H::AffineChartHomotopy{<:Any,N},
    x::AbstractVector,
    t,
) where {N}
    evaluate_and_jacobian!(u, U, H.homotopy, x, t)
    evaluate_chart!(u, H.chart, x)
    jacobian_chart!(U, H.chart, x)
    nothing
end

function taylor!(u, v::Val{N}, H::AffineChartHomotopy, tx::TaylorVector{N}, t) where {N}
    u .= zero(eltype(u))
    taylor!(u, v, H.homotopy, tx, t)
    # affine chart part is always zero since it is an affine linear form
    u
end
