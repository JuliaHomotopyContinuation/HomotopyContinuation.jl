__precompile__()

module PredictorCorrector
    using ..HomConBase
    using ..Homotopy
    import ..HomConBase: HomConAlgorithm, AbstractHomotopy, correct!, predict, solve

    abstract type AbstractPredictorCorrectorHomConAlgorithm <: HomConAlgorithm end

    include("solve.jl")
    include("spherical.jl")
    include("affine.jl")

    export Spherical, Affine
end