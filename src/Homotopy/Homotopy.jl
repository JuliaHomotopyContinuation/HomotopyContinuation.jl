module Homotopy

    using ..HomConBase
    import ..HomConBase: AbstractHomotopy, AbstractPolySystem, evaluate, substitute,
        startsystem, targetsystem,
        differentiate, dt,
        homogenize, ishomogenous, homogenized,
        degrees, weylnorm, nvariables

    import Base: length
    include("straight_line.jl")
    include("gamma_trick.jl")
    include("fixedhomotopy.jl")

    export evaluate, substitute, startsystem, targetsystem, differentiate,
        dt, homogenize, homogenized, degrees, weylnorm, nvariables
    export StraightLineHomotopy, GammaTrickHomotopy, FixedHomotopy
end
