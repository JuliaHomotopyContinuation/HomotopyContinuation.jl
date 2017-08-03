module HomotopyContinuation
    using Reexport

    include("MPoly/MPoly.jl")
    
    export MPoly

    include("HomConBase/HomConBase.jl")

    @reexport using .HomConBase


    include("Homotopy/Homotopy.jl")

    @reexport using .Homotopy

    include("PredictorCorrector/PredictorCorrector.jl")

    @reexport using .PredictorCorrector
end # module
 