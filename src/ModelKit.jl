@reexport module ModelKit

export @var,
    @unique_var,
    Variable,
    Expression,
    coefficients,
    degree,
    degrees,
    differentiate,
    dense_poly,
    evaluate,
    expand,
    exponents_coefficients,
    expressions,
    horner,
    is_homogeneous,
    multi_degrees,
    parameters,
    nparameters,
    nvariables,
    monomials,
    optimize,
    parameters,
    subs,
    support_coefficients,
    rand_poly,
    to_dict,
    to_number,
    variables,
    System,
    Homotopy,
    TaylorVector,
    variable_groups,
    vectors

import LinearAlgebra

include("model_kit/symengine.jl")
include("model_kit/symbolic.jl")
include("model_kit/instructions.jl")
include("model_kit/codegen.jl")

end # module
