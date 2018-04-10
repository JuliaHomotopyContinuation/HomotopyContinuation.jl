using HomotopyContinuation
using Compat.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials
using TestSystems

include("systems_test.jl")
include("new_homotopies_test.jl")
include("utilities_test.jl")
include("problem_test.jl")
include("predictors_test.jl")
include("correctors_test.jl")
include("affine_patches_test.jl")
include("projective_path_tracker_test.jl")
include("pathtracker_iterator_test.jl")
