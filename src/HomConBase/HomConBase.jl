__precompile__()

module HomConBase

    import Base: eltype, gradient


    export AbstractHomotopy
    abstract type AbstractHomotopy{T<:Number} end

    abstract type HomConAlgorithm end
    
    export HomConAlgorithm

    # abstract type AbstractPredictorAlgorithm end
    # abstract type AbstractCorrectorAlgorithm end
    # export AbstractPredictorAlgorithm, AbstractCorrectorAlgorithm
    
    # Predictor Corrector
    function correct! end
    function predict end


    function init end
    function solve end

    
    export init, solve

    include("poly.jl")
    include("polysystem.jl")
    include("homotopy.jl")

    export Poly, evaluate, gradient, homogenize, is_homogenous, deg, nvars, weyl_dot, weyl_norm
    
    export PolySystem, evaluate, jacobian, is_homogenous, homogenize, degrees, nvars, nequations

    export evaluate, startsystem, targetsystem, jacobian, dt, ∂t, nvars, degrees, nequations, is_homogenous, homogenize  # nvars, nequations, degrees, 

end