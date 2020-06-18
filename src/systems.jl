fixed(F::AbstractSystem; kwargs...) = F
fixed(F::System; compile::Bool = COMPILE_DEFAULT[], kwargs...) =
    compile ? CompiledSystem(F; kwargs...) : InterpretedSystem(F; kwargs...)

include("systems/fixed_parameter_system.jl")
include("systems/randomized_system.jl")
include("systems/affine_chart_system.jl")
include("systems/composition_system.jl")
include("systems/start_pair_system.jl")
