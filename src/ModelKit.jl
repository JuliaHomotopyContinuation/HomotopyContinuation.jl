@reexport module ModelKit

export @var,
    @unique_var,
    AbstractHomotopy,
    AbstractSystem,
    CompiledHomotopy,
    CompiledSystem,
    Expression,
    Homotopy,
    InterpretedHomotopy,
    InterpretedSystem,
    System,
    TaylorVector,
    Variable,
    coefficients,
    degree,
    degrees,
    differentiate,
    dense_poly,
    evaluate,
    evaluate!,
    evaluate_and_jacobian!,
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
    taylor!,
    to_dict,
    to_number,
    variables,
    variable_groups,
    vectors

using ..DoubleDouble: ComplexDF64
import LinearAlgebra
import MultivariatePolynomials:
    MultivariatePolynomials,
    coefficients,
    degree,
    differentiate,
    nvariables,
    monomials,
    subs,
    variables
using Parameters: @unpack
const MP = MultivariatePolynomials

include("model_kit/symengine.jl")
include("model_kit/symbolic.jl")
include("model_kit/taylor_vector.jl")
include("model_kit/instructions.jl")
include("model_kit/abstract_types.jl")
include("model_kit/compiler.jl")
include("model_kit/compiled.jl")
include("model_kit/interpreter.jl")
include("model_kit/interpreted.jl")

end # module
