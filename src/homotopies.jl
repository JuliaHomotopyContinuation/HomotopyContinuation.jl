export AbstractHomotopy

"""
    AbstractHomotopy

An abstract type representing a homotopy ``H(x,t)``.

The following homotopies are available:

* [`AffineChartHomotopy`](@ref)
* [`AffineSubspaceHomotopy`](@ref)
* [`ModelKitHomotopy`](@ref)
* [`ParameterHomotopy`](@ref)
* [`StraightLineHomotopy`](@ref)
"""
abstract type AbstractHomotopy end

Base.size(H::AbstractHomotopy, i::Integer) = size(H)[i]
function set_solution!(x::AbstractVector, H::AbstractHomotopy, y::AbstractVector, t)
    x .= y
end
get_solution(H::AbstractHomotopy, x::AbstractVector, t) = copy(x)

# internal only
include("homotopies/toric_homotopy.jl")

# public, these should be linked on the top
include("homotopies/affine_chart_homotopy.jl")
include("homotopies/affine_subspace_homotopy.jl")
include("homotopies/model_kit_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")


include("homotopies/numerical_differentiation.jl")

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
    ND::NumericalDifferentiation,
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
