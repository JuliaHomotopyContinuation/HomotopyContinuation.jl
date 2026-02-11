using HomotopyContinuation
using LinearAlgebra, Test, Parameters, Random
using Arblib
using HomotopyContinuation.DoubleDouble: ComplexDF64
const HC = HomotopyContinuation

set_default_compile(:none)

include("test_systems.jl")

Random.seed!(0x8b868a97)

include("./model_kit/symbolic_test.jl")
include("./model_kit/operations_test.jl")
include("./model_kit/e2e_test.jl")
include("./model_kit/slp_test.jl")

include("systems_test.jl")
include("tracker_test.jl")
include("endgame_tracker_test.jl")
include("polyhedral_test.jl")
include("solve_test.jl")
include("result_test.jl")
include("endgame_test.jl")
