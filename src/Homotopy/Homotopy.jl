module Homotopy

    using ..HomConBase
    import ..HomConBase: AbstractHomotopy, evaluate,
        startsystem, targetsystem,
        differentiate, dt,
        homogenize, ishomogenous, homogenized,
        degrees, weylnorm, nvariables

    import Base: length
    include("straight_line.jl")
    include("gamma_trick.jl")

    export evaluate, startsystem, targetsystem, differentiate,
        dt, homogenize, homogenized, degrees, weylnorm, nvars, nequations
    export StraightLineHomotopy, GammaTrickHomotopy
end
