__precompile__()

module PredictorCorrector
    using ..HomConBase
    using ..Homotopy
    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    import ..HomConBase: HomConAlgorithm, AbstractHomotopy, correct!, predict, solve

    abstract type AbstractPredictorCorrectorHomConAlgorithm <: HomConAlgorithm end
    
    export AbstractPredictorCorrectorHomConAlgorithm

    # These functions are used in `solve`, you want to overwrite if necessary
    prepare_homotopy(H::AbstractHomotopy, ::AbstractPredictorCorrectorHomConAlgorithm) = H
    prepare_start_value(start_value, ::AbstractHomotopy, ::AbstractPredictorCorrectorHomConAlgorithm) = start_value
    hom_var_index(::AbstractHomotopy, ::AbstractPredictorCorrectorHomConAlgorithm) = Nullable{Int64}()

    export prepare_homotopy, prepare_start_value, hom_var_index

    include("solve.jl")
    include("spherical.jl")
    include("affine.jl")

    export Spherical, Affine
end