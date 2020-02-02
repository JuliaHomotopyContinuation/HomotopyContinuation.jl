module ModelKit

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
    expand,
    System,
    Homotopy,
    evaluate!,
    jacobian!,
    evaluate_and_jacobian!,
    diff_t!,
    evaluate,
    jacobian,
    evaluate_and_jacobian,
    diff_t

import LinearAlgebra

include("model_kit/symengine.jl")
include("model_kit/symbolic.jl")
include("model_kit/instructions.jl")
include("model_kit/codegen.jl")

end # module
