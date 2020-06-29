"""
    fixed(H::Homotopy; compile::Bool = COMPILE_DEFAULT[])

Constructs either a [`CompiledHomotopy`](@ref) (if `compile = true`) or an
[`InterpretedHomotopy`](@ref) (if `compile = false`).
"""
fixed(H::Homotopy; compile::Bool = COMPILE_DEFAULT[], kwargs...) =
    compile ? CompiledHomotopy(H; kwargs...) : InterpretedHomotopy(H; kwargs...)
fixed(H::AbstractHomotopy; kwargs...) = H


function set_solution!(x::AbstractVector, H::AbstractHomotopy, y::AbstractVector, t)
    x .= y
end
get_solution(H::AbstractHomotopy, x::AbstractVector, t) = copy(x)

start_parameters!(H::AbstractHomotopy, p) = H
target_parameters!(H::AbstractHomotopy, p) = H
# internal only
include("homotopies/toric_homotopy.jl")

# public, these should be linked on the top
include("homotopies/affine_chart_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/coefficient_homotopy.jl")
include("homotopies/extrinsic_subspace_homotopy.jl")
include("homotopies/intrinsic_subspace_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")
include("homotopies/fixed_parameter_homotopy.jl")
