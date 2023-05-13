export AffineChartHomotopy, on_affine_chart

"""
    AffineChartHomotopy(H::AbstractHomotopy, v::PVector{T,N})

Given a homotopy ``H(x,t): (ℙ^{m_1} × ⋯ × ℙ^{m_N}) × ℂ → ℂⁿ`` this creates a new affine
homotopy ``H̄`` which operates on the affine chart defined by the vector
``v ∈ ℙ^{m_1} × ⋯ × ℙ^{m_N}`` and the augmented conditions ``vᵀx = 1``.
"""
# struct AffineChartHomotopy{H<:AbstractHomotopy,N} <: AbstractHomotopy
#     homotopy::H
#     chart::PVector{ComplexF64,N}
# end
struct AffineChartHomotopy{T<:AbstractHomotopy} <: AbstractHomotopy
    H::T
    composed::SystemHomotopy{AffineChartSystem{HomotopySystem{T}}}
end
AffineChartHomotopy(
    H::AbstractHomotopy,
    chart::PVector;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = AffineChartHomotopy(
    H,
    system_as_homotopy(AffineChartSytem(homotopy_as_system(H), chart)),
)

Base.size(H::AffineChartHomotopy) = size(H.composed)

(H::AffineChartHomotopy)(x, t, p::Nothing = nothing) = H.composed(x, t)

function set_solution!(x::AbstractVector, H::AffineChartHomotopy, y::AbstractVector, t)
    set_solution!(x, H.composed, y, t)
end

start_parameters!(H::AffineChartHomotopy, p) = start_parameters!(H.composed, p)
target_parameters!(H::AffineChartHomotopy, p) = target_parameters!(H.composed, p)
parameters!(H::AffineChartHomotopy, p, q) = parameters!(H.composed, p, q)

function ModelKit.evaluate!(u, H::AffineChartHomotopy, x, t, p::Nothing = nothing)
    evaluate!(u, H.composed, x, t)
end
function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::AffineChartHomotopy,
    x,
    t,
    p::Nothing = nothing,
)
    evaluate_and_jacobian!(u, U, H.composed, x, t)
end

function ModelKit.taylor!(u, v::Val, H::AffineChartHomotopy, tx, tt, p::Nothing = nothing)
    taylor!(u, v, H.composed, tx, tt)
end



