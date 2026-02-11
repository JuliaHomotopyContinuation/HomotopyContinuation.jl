using HomotopyContinuation
using LinearAlgebra, Test, Parameters, Random
using HomotopyContinuation.DoubleDouble: ComplexDF64
using ProjectiveVectors: PVector, affine_chart
const HC = HomotopyContinuation

set_default_compile(:none)

include("test_systems.jl")

hc_testset = get(ENV, "HC_TESTSET", "none")

Random.seed!(0x8b868a97)

@testset "HomotopyContinuation" begin
    if (hc_testset == "none") || (hc_testset == "model_kit")
        # model_kit tests are intentionally excluded in this minimal build
    end
    if (hc_testset == "none") || (hc_testset == "non_model_kit")
        include("systems_test.jl")
        include("tracker_test.jl")
        include("endgame_tracker_test.jl")
        include("polyhedral_test.jl")
        include("solve_test.jl")
        include("result_test.jl")
        include("endgame_test.jl")
    end

end
