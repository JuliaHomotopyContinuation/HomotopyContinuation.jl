using HomotopyContinuation
using Compat.Test
using TestSystems

# We order the tests such that isolated things are tested first
include("utilities_test.jl")
include("systems_test.jl")
include("homotopies_test.jl")
include("problem_test.jl")
include("predictors_test.jl")
include("correctors_test.jl")
include("affine_patches_test.jl")
include("projective_path_tracker_test.jl")
include("pathtracker_iterator_test.jl")
include("solve_test.jl")
include("integration_tests.jl")
