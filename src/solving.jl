module Solving

using Compat
using LinearAlgebra

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ProgressMeter
# import Juno

import ..AffinePatches
import ..Endgaming
import ..Homotopies
import ..PathTracking
import ..Parallel
import ..Predictors
import ..Problems
import ..PatchSwitching
import ..StepLength
import ..Systems

export Solver,
    solve,
    multiplicities,
    Result

include("solving/multiplicities.jl")
include("solving/path_result.jl")
include("solving/result.jl")
include("solving/types.jl")
include("solving/solve.jl")


end
