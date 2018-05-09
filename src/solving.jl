module Solving

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ProgressMeter
import Juno

import ..Endgame
import ..Homotopies
import ..PathTracking
import ..Parallel
import ..Problems
import ..Systems

export Solver,
    solve,
    Result,
    real_solutions

include("solving/path_result.jl")
include("solving/result.jl")
include("solving/types.jl")
include("solving/solve.jl")


end
