using HomotopyContinuation2, LinearAlgebra, Test
const HC2 = HomotopyContinuation2

@testset "HomotopyContinuation2" begin
    include("utils_test.jl")
    include("double_double_test.jl")
    include("norm_test.jl")
    include("model_kit_test.jl")
    include("linear_algebra_test.jl")
    include("homotopies_test.jl")
    include("tracker_test.jl")
end
