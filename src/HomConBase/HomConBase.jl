module HomConBase

    #using Reexport
    import Base: eltype, gradient

    import FixedPolySystem: Poly, PolySystem,
        nvariables, variables, polynomials, degrees,
        evaluate, evaluate!, differentiate,
        ishomogenous, homogenize, homogenized, dehomogenize,
        weyldot, weylnorm
    # Exports from FixedPolySystem
    export Poly, PolySystem,
        nvariables, variables, polynomials, degrees,
        evaluate, evaluate!, differentiate,
        ishomogenous, homogenize, homogenized, dehomogenize,
        weyldot, weylnorm


    export AbstractHomotopy
    abstract type AbstractHomotopy{T<:Number} end

    abstract type HomConAlgorithm end

    export HomConAlgorithm


    function init end
    function solve end


    export init, solve

    include("homotopy.jl")
    include("utilities.jl")

    export evaluate, startsystem, targetsystem,
        differentiate, dt, ∂t, nvariables, degrees,
        nequations, ishomogenous, homogenize, homogenized, weylnorm

    export affine, projective, totaldegree
end
