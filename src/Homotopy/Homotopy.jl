__precompile__()

module Homotopy

    using ..HomConBase
    import ..HomConBase: AbstractHomotopy, evaluate, startsystem, targetsystem, jacobian,
        dt, homogenize, degrees, weyl_norm, nvars, nequations
    include("straight_line.jl")
    include("gamma_trick.jl")

    export StraightLineHomotopy, GammaTrickHomotopy
end