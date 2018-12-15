module Problems

import DynamicPolynomials
import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import Random

import ..Input
import ..Input: AbstractInput, TotalDegree, StartTarget, ParameterSystem, MPPolys
import ..Homotopies: AbstractHomotopy, StraightLineHomotopy, ParameterHomotopy
import ProjectiveVectors
import ProjectiveVectors: PVector
import ..Systems: AbstractSystem, SPSystem, FPSystem
import ..Systems

using ..Utilities

include("problems/homogenization.jl")
include("problems/problem.jl")
include("problems/totaldegree_solutions.jl")
include("problems/problem_startsolutions.jl")

end
