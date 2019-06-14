module HomotopyContinuation

    import DynamicPolynomials
    import ElasticArrays
    import FixedPolynomials
    import LinearAlgebra
    import MixedSubdivisions
    import MultivariatePolynomials
    import PrettyTables
    import Printf
    import ProgressMeter
    import ProjectiveVectors
    import Random
    import StaticArrays
    import StaticPolynomials
    import TreeViews

    import Parameters: @pack!, @unpack
    import DynamicPolynomials: @polyvar
    import ProjectiveVectors: PVector
    import StaticArrays: SVector, @SVector
    import Test: @test

    const FP = FixedPolynomials
    const MP = MultivariatePolynomials
    const SP = StaticPolynomials

    export @polyvar

    include("utilities.jl")
    include("affine_patches.jl")

    include("systems_and_homotopies.jl")
    include("input.jl")
    include("problems.jl")
    include("totaldegree.jl")

    include("newton.jl")
    include("predictors.jl")
    include("correctors.jl")

    include("core_tracker.jl")
    include("path_tracker.jl")
    include("polyhedral.jl")
    include("solve.jl")
    include("monodromy.jl")

    import LinearAlgebra: issuccess
    export issuccess
end
