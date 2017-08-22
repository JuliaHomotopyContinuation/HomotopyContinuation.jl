module HomConBase

    #using Reexport
    import Base: eltype

    import FixedPolySystem: AbstractPolySystem, Poly, PolySystem,
        nvariables, variables, removepoly, degrees,
        evaluate, evaluate!, substitute, differentiate,
        ishomogenous, homogenize, homogenized, dehomogenize,
        weyldot, weylnorm, coefftype
    # Exports from FixedPolySystem
    export AbstractPolySystem, Poly, PolySystem,
        nvariables, variables, polynomials, degrees,
        evaluate, evaluate!, differentiate, substitute, removepoly,
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

    export evaluate, substitute, startsystem, targetsystem,
        differentiate, dt, ∂t, nvariables, degrees,
        nequations, removepoly, ishomogenous, homogenize, homogenized, weylnorm, coefftype

    export affine, projective, totaldegree
end
