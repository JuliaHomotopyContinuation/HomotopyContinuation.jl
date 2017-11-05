using HomotopyContinuation
using Base.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials

include("pathtracker_test.jl")
include("affine_test.jl")
include("solver_test.jl")
include("spherical_test.jl")
include("cauchyendgame_test.jl")
#include("result_test.jl")
include("test_systems_test.jl")
include("solve_test.jl")
include("patch_test.jl")
