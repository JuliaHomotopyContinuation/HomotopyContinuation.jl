module HomotopyContinuation

    import DynamicPolynomials
    import LinearAlgebra
    import MultivariatePolynomials
    import ProgressMeter
    import ProjectiveVectors
    import Random
    import TreeViews

    import DynamicPolynomials: @polyvar
    import ProjectiveVectors: PVector
    import Test: @test

    const MP = MultivariatePolynomials

    export @polyvar
    export AffinePatches,
        Utilities

    include("utilities.jl")
    include("parallel.jl")
    include("affine_patches.jl")

    include("systems_and_homotopies.jl")
    include("input.jl")
    include("problems.jl")
    include("totaldegree.jl")
    include("predictors.jl")
    include("correctors.jl")

    include("path_tracking.jl")
    include("endgaming.jl")

    include("solving.jl")
    include("solve.jl")
    include("monodromy.jl")

    import .Utilities: homogenize, uniquevar, ishomogenous
    export homogenize, uniquevar, ishomogenous

    import LinearAlgebra: issuccess
    export issuccess
end
