export fixed

"""
    fixed(F::System; compile::Bool = COMPILE_DEFAULT[])

Constructs either a [`CompiledSystem`](@ref) (if `compile = true`) or an
[`InterpretedSystem`](@ref) (if `compile = false`).
"""
fixed(F::System; compile::Bool = COMPILE_DEFAULT[], kwargs...) =
    compile ? CompiledSystem(F; kwargs...) : InterpretedSystem(F; kwargs...)
fixed(F::AbstractSystem; kwargs...) = F


include("systems/fixed_parameter_system.jl")
include("systems/randomized_system.jl")
include("systems/affine_chart_system.jl")
include("systems/composition_system.jl")
include("systems/start_pair_system.jl")
