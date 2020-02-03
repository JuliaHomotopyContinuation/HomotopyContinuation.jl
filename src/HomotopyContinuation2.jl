module HomotopyContinuation2

export ModelKit

import LinearAlgebra
import StructArrays

const LA = LinearAlgebra

include("utils.jl")

include("ModelKit.jl")
import .ModelKit

using .ModelKit: @var,
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
                 Homotopy

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
       Homotopy

include("linear_algebra.jl")
include("homotopies.jl")

end
