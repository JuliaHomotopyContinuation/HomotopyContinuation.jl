@reexport module ModelKit

export @var,
    @unique_var,
    Variable,
    Expression,
    variables,
    nvariables,
    subs,
    evaluate,
    differentiate,
    monomials,
    rand_poly,
    coefficients,
    exponents_coefficients,
    expand,
    System,
    Homotopy

import LinearAlgebra

include("model_kit/symengine.jl")
include("model_kit/symbolic.jl")
include("model_kit/instructions.jl")
include("model_kit/codegen.jl")

end # module
