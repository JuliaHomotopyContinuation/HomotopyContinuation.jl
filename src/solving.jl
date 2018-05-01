module Solving

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ..Endgame
import ..Homotopies
import ..PathTracking
import ..Parallel
import ..Problems
import ..Systems

export Solver,
    solve


include("solving/path_result.jl")
include("solving/types.jl")
include("solving/solve.jl")

end
