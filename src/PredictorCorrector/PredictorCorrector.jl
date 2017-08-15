module PredictorCorrector
    using ..HomConBase
    using ..Homotopy
    import ..HomConBase: HomConAlgorithm, AbstractHomotopy, solve

    abstract type AbstractPredictorCorrectorHomConAlgorithm <: HomConAlgorithm end

    export AbstractPredictorCorrectorHomConAlgorithm

    function correct! end
    function predict end

    export correct!, predict

    include("solve.jl")
    include("spherical.jl")
    include("affine.jl")

    export Spherical, Affine, solve
end
