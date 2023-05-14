# internal only
include("homotopies/t_embedding.jl")
include("homotopies/general_parameter_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")
include("homotopies/system_homotopy_stack.jl")
include("homotopies/extrinsic_linear_subspace_homotopy.jl")
include("homotopies/system_homotopy.jl")
include("homotopies/system_homotopy_composition.jl")
include("homotopies/instrinsic_linear_subspace_homotopy.jl")
include("homotopies/linear_subspace_geodesic_homotopy.jl")
include("homotopies/mixed_homotopy.jl")
include("homotopies/affine_chart_homotopy.jl")


# include("homotopies/toric_homotopy.jl")


"""
    fixed(H::Homotopy; compile::Union{Bool,Symbol} = $(COMPILE_DEFAULT[]))

Constructs either a [`CompiledHomotopy`](@ref) (if `compile = :all`), an
[`InterpretedHomotopy`](@ref) (if `compile = :none`) or a
[`MixedHomotopy`](@ref) (`compile = :mixed`).
"""
function fixed(H::Homotopy; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[], kwargs...)
    if compile == true || compile == :all
        CompiledHomotopy(H; kwargs...)
    elseif compile == false || compile == :none
        InterpretedHomotopy(H; kwargs...)
    elseif compile == :mixed
        MixedHomotopy(H; kwargs...)
    else
        error("Unknown argument $compile for keyword `compile`.")
    end
end
ModelKit.variable_groups(::AbstractHomotopy) = nothing
fixed(H::AbstractHomotopy; kwargs...) = H

function set_solution!(x::AbstractVector, H::AbstractHomotopy, y::AbstractVector, t)
    x .= y
end
get_solution(H::AbstractHomotopy, x::AbstractVector, t) = copy(x)

start_parameters!(H::AbstractHomotopy, p) = H
target_parameters!(H::AbstractHomotopy, p) = H

ModelKit.taylor!(u, v::Val, H::AbstractHomotopy, tx, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)
