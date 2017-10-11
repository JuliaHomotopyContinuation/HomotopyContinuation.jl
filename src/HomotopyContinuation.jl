__precompile__()

module HomotopyContinuation

    using Reexport
    @reexport using Homotopy

    import MultivariatePolynomials
    const MP = MultivariatePolynomials

    abstract type AbstractHomotopyContinuationAlgorithm end
    export AbstractHomotopyContinuationAlgorithm

    include("predictorcorrector.jl")
    include("solve.jl")
    include("testsystems.jl")

    export TestSystems
end # module
