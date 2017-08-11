__precompile__()

module Homotopy

    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    using ..HomConBase
    import ..HomConBase: AbstractHomotopy, evaluate, startsystem, targetsystem, jacobian,
        dt, homogenize, degrees, weyl_norm, nvars, vars, nequations
    include("straight_line.jl")
    include("gamma_trick.jl")

    export StraightLineHomotopy, GammaTrickHomotopy
end