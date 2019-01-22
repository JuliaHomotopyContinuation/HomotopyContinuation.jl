using LinearAlgebra

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ProgressMeter
import TreeViews

import ..AffinePatches
import ..Homotopies
import ..Parallel
import ..Predictors
import ..Problems
import ..Systems

using ..Utilities

export Solver,
    solve,
    multiplicities,
    Result

include("solving/multiplicities.jl")
include("solving/path_result.jl")
include("solving/result.jl")
include("solving/solver.jl")
include("solving/solve.jl")
