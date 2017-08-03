module PredictorCorrector
    using ..HomConBase
    using ..Homotopy
    import ..HomConBase: HomConAlgorithm, AbstractHomotopy, correct!, predict, solve

    abstract type AbstractPredictorCorrectorHomConAlgorithm <: HomConAlgorithm end

    include("solve.jl")
    include("spherical.jl")
end