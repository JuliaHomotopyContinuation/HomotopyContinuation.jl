module Homotopy

    using ..HomConBase
    import ..HomConBase: AbstractHomotopy, AbstractPolySystem, evaluate,
        startsystem, targetsystem,
        differentiate, dt,
        homogenize, ishomogenous, homogenized,
        degrees, weylnorm, nvariables

    import Base: length
    include("straight_line.jl")
    include("gamma_trick.jl")
    include("fixedhomotopy.jl")

    export evaluate, startsystem, targetsystem, differentiate,
        dt, homogenize, homogenized, degrees, weylnorm, nvars, nequations
    export StraightLineHomotopy, GammaTrickHomotopy, FixedHomotopy
end
