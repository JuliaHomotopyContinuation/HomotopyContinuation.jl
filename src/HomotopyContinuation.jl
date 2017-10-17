__precompile__()

module HomotopyContinuation

    using Reexport
    @reexport using Homotopy

    using Parameters

    using HigherPrecision

    import MultivariatePolynomials
    const MP = MultivariatePolynomials


    abstract type AbstractPathtrackingAlgorithm end
    abstract type AbstractEndgameAlgorithm end

    abstract type AbstractPathtrackerCache{T<:Complex} end
    abstract type AbstractEndgameCache end

    export AbstractPathtrackingAlgorithm, AbstractEndgameAlgorithm,
        AbstractPathtrackerCache, AbstractEndgameCache

    include("patches.jl")
    include("utilities.jl")

    include("pathtracker/type.jl")
    include("pathtracker/initialize_modify.jl")
    include("pathtracker/interface.jl")
    include("pathtracker/result.jl")

    include("pathtracking_algorithms/spherical.jl")
    include("pathtracking_caches/spherical_cache.jl")

    include("endgamer/type.jl")
    include("endgamer/initialize_modify.jl")
    include("endgamer/iterator_interface.jl")
    include("endgamer/result.jl")
    include("endgame_algorithms/cauchy.jl")
    include("endgame_cache/cauchy_cache.jl")

    include("solver.jl")

    include("result.jl")
    # include("predictorcorrector.jl")

    include("solve.jl")
    include("testsystems.jl")

    export TestSystems
end # module
