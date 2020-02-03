module HomotopyContinuation2

export ModelKit

using DoubleFloats: Double64
import LinearAlgebra
using Parameters: @unpack
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

include("norm.jl")
include("linear_algebra.jl")
include("homotopies.jl")
include("predictors.jl")
include("newton_corrector.jl")
include("tracker.jl")

end
