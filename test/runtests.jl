using HomotopyContinuation2, LinearAlgebra, Test, StaticArrays, Parameters
const HC2 = HomotopyContinuation2

include("test_systems.jl")

@testset "HomotopyContinuation2" begin
    include("utils_test.jl")
    include("double_double_test.jl")
    include("binomial_system_test.jl")
    include("norm_test.jl")
    include("model_kit_test.jl")
    include("linear_algebra_test.jl")
    include("homotopies_test.jl")
    include("tracker_test.jl")
    include("valuation_test.jl")
    include("path_tracker_test.jl")
    include("endgame_test.jl")
    include("polyhedral_test.jl")
    include("linear_test.jl")

    if "extensive" in ARGS
        include("extensive/extensive_test.jl")
    end
end
