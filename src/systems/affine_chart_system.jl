export AffineChartSystem, on_affine_chart

"""
    AffineChartSystem(F::AbstractSystem, v::PVector{T,N})

Given a system ``F(x): (ℙ^{m_1} × ⋯ × ℙ^{m_N}) → ℂⁿ`` this creates a new affine
system ``F′`` which operates on the affine chart defined by the vector
``v ∈ ℙ^{m_1} × ⋯ × ℙ^{m_N}`` and the augmented conditions ``vᵀx = 1``.
"""
struct AffineChartSystem{S<:AbstractSystem,N} <: AbstractSystem
    system::S
    chart::PVector{ComplexF64,N}
    stacked::StackedSystem{S,LinearSystem}
end

function AffineChartSystem(system::AbstractSystem, v::PVector{ComplexF64,N}) where {N}
    m, n = size(system)
    A = zeros(ComplexF64, N, n)
    b = ones(ComplexF64, N)
    ranges = dimension_indices(v)
    for (k, range) in enumerate(ranges)
        for i in range
            A[k, i] = v[i]
        end
    end
    L = LinearSystem(A, b; variables = variables(system))
    stacked = StackedSystem(system, L)

    AffineChartSystem(system, v, stacked)
end

ModelKit.variables(F::AffineChartSystem) = variables(F.system)
ModelKit.parameters(F::AffineChartSystem) = parameters(F.system)
ModelKit.variable_groups(F::AffineChartSystem) = variable_groups(F.system)

"""
    on_affine_chart(F::Union{System,AbstractSystem}, dimensions)

Construct an `AffineChartSystem` on a randomly generated chart `v`. Each entry is drawn
idepdently from a univariate normal distribution.
"""
on_affine_chart(
    F::System,
    dims = nothing;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...,
) = on_affine_chart(fixed(F; compile = compile), dims)
function on_affine_chart(
    F::AbstractSystem,
    dims = nothing,
    compile::Union{Bool,Symbol} = true,
)
    vargroups = variable_groups(F)
    if vargroups === nothing
        dims = [size(F, 2) - 1]
    else
        dims = length.(vargroups,) .- 1
    end
    chart = PVector(randn(ComplexF64, sum(dims) + length(dims)), tuple(dims...))
    AffineChartSystem(F, chart)
end

function Base.size(F::AffineChartSystem{<:AbstractSystem,N}) where {N}
    m, n = size(F.system)
    (m + N, n)
end

function set_solution!(x, F::AffineChartSystem, y)
    x .= y
    on_chart!(x, F.chart)
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

(F::AffineChartSystem)(x, p = nothing) = F.stacked(x, p)
ModelKit.evaluate!(u, F::AffineChartSystem, x::AbstractVector, p = nothing) =
    evaluate!(u, F.stacked, x, p)
ModelKit.evaluate_and_jacobian!(u, U, F::AffineChartSystem, x, p = nothing) =
    evaluate_and_jacobian!(u, U, F.stacked, x, p)
ModelKit.taylor!(u, v::Val, F::AffineChartSystem, tx, p = nothing) =
    taylor!(u, v, F.stacked, tx, p)
