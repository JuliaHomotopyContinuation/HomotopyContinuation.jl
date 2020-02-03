using HomotopyContinuation2, LinearAlgebra, Test
const HC2 = HomotopyContinuation2

@testset "HomotopyContinuation2" begin
    include("model_kit_test.jl")
    include("linear_algebra_test.jl")
    include("homotopies_test.jl")
end
