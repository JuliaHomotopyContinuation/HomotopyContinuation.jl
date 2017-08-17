using HomotopyContinuation
using Base.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials

include("base_utilities_test.jl")

include("straight_line_test.jl")
include("gamma_trick_test.jl")
include("spherical_test.jl")
include("affine_test.jl")
include("test_systems_test.jl")
