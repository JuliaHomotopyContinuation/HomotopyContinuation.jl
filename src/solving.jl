module Solving

using Compat

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ProgressMeter
import Juno

import ..AffinePatches
import ..Endgame
import ..Homotopies
import ..PathTracking
import ..Parallel
import ..Problems
import ..Systems

export Solver,
    solve,
    Result,
    finite_solutions,
    solutions_at_infinity,
    singular_solutions,
    failed_paths,
    real_solutions

include("solving/path_result.jl")
include("solving/result.jl")
include("solving/types.jl")
include("solving/solve.jl")


end
