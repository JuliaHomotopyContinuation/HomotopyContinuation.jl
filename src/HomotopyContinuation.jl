module HomotopyContinuation
    include("HomConBase/HomConBase.jl")

    using .HomConBase

    # Exports defined in HomConBase
    export AbstractHomotopy
    export HomConAlgorithm
    export init, solve
    export evaluate, startsystem, targetsystem,
        differentiate, dt, ∂t, nvariables, degrees,
        nequations, ishomogenous, homogenize
    export affine, projective, totaldegree

    # Exports from FixedPolySystem in HomConBase
    export Poly, PolySystem,
        nvariables, variables, polynomials,
        evaluate, evaluate!, differentiate,
        ishomogenous, homogenize, homogenized, dehomogenize,
        weyldot, weylnorm

    include("Homotopy/Homotopy.jl")

    using .Homotopy
    # Homotopy Exports
    export evaluate, startsystem, targetsystem, differentiate,
        dt, homogenize, homogenized, degrees, weylnorm, nvars, nequations
    export StraightLineHomotopy, GammaTrickHomotopy


    include("PredictorCorrector/PredictorCorrector.jl")

    using .PredictorCorrector
    #PredictorCorrector exports
    export PredictorCorrector
    export correct!, predict
    export Spherical, Affine, trackpath

    #
    include("TestSystems/TestSystems.jl")
    export TestSystems
end # module
