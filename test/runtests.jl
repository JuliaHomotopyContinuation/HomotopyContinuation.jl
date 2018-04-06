#addprocs(4)
using HomotopyContinuation
using Compat.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials


include("pathtracker_test.jl")
include("solver_test.jl")
include("affine_test.jl")
include("spherical_test.jl")
include("cauchyendgame_test.jl")
include("pathcrossing_test.jl")
#include("result_test.jl")
include("test_systems_test.jl")
include("solve_test.jl")
# include("patch_test.jl")

using TestSystems
include("systems_test.jl")
include("new_homotopies_test.jl")
include("utilities_test.jl")
include("problem_test.jl")
