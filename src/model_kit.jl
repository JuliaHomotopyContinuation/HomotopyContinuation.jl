# These functions should be overloaded from ModelKit, so we need to make sure
# that they are already defined in the HomotopyContinuation module.
function evaluate end
function evaluate! end
function evaluate_and_jacobian end
function evaluate_and_jacobian! end
function jacobian end
function jacobian! end
function dt end
function dt! end
function jacobian_and_dt end
function jacobian_and_dt! end


module ModelKit

using OrderedCollections: OrderedDict
using StaticArrays: @SVector, @SMatrix

import LinearAlgebra: det, dot
import Latexify
# Overload these functions vom MP to be able to export this without any name clashes
import MultivariatePolynomials: variables, differentiate, subs, monomials

export Expression, Constant, Variable, Operation
import ..HomotopyContinuation: evaluate,
                               evaluate!,
                               evaluate_and_jacobian,
                               evaluate_and_jacobian!,
                               jacobian,
                               jacobian!,
                               dt,
                               dt!,
                               jacobian_and_dt,
                               jacobian_and_dt!

export @var,
       @unique_var,
       subs,
       variables,
       differentiate,
       monomials,
       Compiled,
       CompiledSystem,
       CompiledHomotopy,
       compile,
       interpreted,
       System,
       Homotopy,
       # reexport
       evaluate,
       evaluate!,
       evaluate_gradient,
       evaluate_and_jacobian,
       evaluate_and_jacobian!,
       jacobian,
       jacobian!,
       dt,
       dt!,
       jacobian_and_dt,
       jacobian_and_dt!

include("model_kit/expression.jl")
include("model_kit/codegen.jl")

end
