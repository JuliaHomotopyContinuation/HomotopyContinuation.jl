@reexport module ModelKit

export @var,
    @unique_var,
    AbstractSystem,
    AbstractHomotopy,
    Expression,
    Variable,
    Interpreter,
    System,
    InterpretedSystem,
    CompiledSystem,
    Homotopy,
    InterpretedHomotopy,
    CompiledHomotopy,
    TaylorVector,
    TruncatedTaylorSeries,
    Variable,
    coefficients,
    coeffs_as_dense_poly,
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
    interpreter,
    is_homogeneous,
    jacobian,
    jacobian!,
    multi_degrees,
    parameters,
    nparameters,
    nvariables,
    monomials,
    optimize,
    parameters,
    poly_from_coefficients_exponents,
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

import Arblib: Arblib, Acb, AcbRef, AcbRefVector
import SimpleGraphs
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

include("./model_kit/symengine.jl")
include("./model_kit/symbolic.jl")
include("./model_kit/operations.jl")
include("./model_kit/intermediate_representation.jl")
include("./model_kit/taylor.jl")
include("./model_kit/acb.jl")
include("./model_kit/instruction_sequence.jl")
include("./model_kit/instruction_interpreter.jl")
include("./model_kit/abstract_system_homotopy.jl")
include("./model_kit/interpreted_system.jl")
include("./model_kit/interpreted_homotopy.jl")
include("./model_kit/compiled_system_homotopy.jl")

end # module
