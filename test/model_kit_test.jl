using HomotopyContinuation.ModelKit

@testset "ModelKit" begin
    include("./model_kit/symbolic_test.jl")
    include("./model_kit/operations_test.jl")
    include("./model_kit/interpreter_dispatch_test.jl")
    include("./model_kit/interpreted_precision_test.jl")
    include("./model_kit/incremental_taylor_test.jl")
    include("./model_kit/double_double_operations_test.jl")
    include("./model_kit/e2e_test.jl")
    include("./model_kit/slp_test.jl")
end
