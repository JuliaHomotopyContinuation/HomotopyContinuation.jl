module PredictorCorrector
    using ..HomConBase
    using ..Homotopy
    import ..HomConBase: HomConAlgorithm, AbstractHomotopy, solve

    export AbstractPredictorCorrectorAlgorithm
    abstract type AbstractPredictorCorrectorAlgorithm{Projective} <: HomConAlgorithm end

    const APCA{P} = AbstractPredictorCorrectorAlgorithm{P}

    """
        predict(alg, H, J_H, ∂H∂t, x, t, Δt)

    Make a prediction step for the algorithm `alg` and the homotopy `H` with the jacobian
    `J_H`, the derivative w.r.t t `∂H∂t` at `x` to the time `t` with a stepwidth `Δt`.
    """
    function predict end

    """
        correct!(u, alg, H, J_H, x, t, tol, max_iterations)::Bool

    Make a correction step for the algorithm `alg` and the homotopy `H` with the jacobian
    `J_H` at `x` to the time `t`. Stores the result in `u` and returns a boolean indicating
    whether ``|u-x|<tol`` withing `max_iterations` iterations.
    """
    function correct! end

    export correct!, predict



    include("spherical.jl")
    include("affine.jl")

    include("prepare.jl")
    include("trackpath.jl")
    include("result_types.jl")
    include("cauchyendgame.jl")
    include("solve.jl")

    export Spherical, Affine, solve, trackpath
end
