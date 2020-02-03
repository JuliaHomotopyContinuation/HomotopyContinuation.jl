abstract type AbstractHomotopy end

Base.size(H::AbstractHomotopy, i::Integer) = size(H)[i]

include("homotopies/model_kit_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
