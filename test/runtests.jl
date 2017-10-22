using HomotopyContinuation
using Base.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials

#include("affine_test.jl")
include("spherical_test.jl")
include("cauchyendgame_test.jl")
#include("result_test.jl")
include("test_systems_test.jl")
include("solve_test.jl")
