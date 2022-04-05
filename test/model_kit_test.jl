using HomotopyContinuation.ModelKit

@testset "ModelKit" begin
    include("./model_kit/symbolic_test.jl")
    include("./model_kit/operations_test.jl")
    include("./model_kit/e2e_test.jl")
    include("./model_kit/slp_test.jl")
end
