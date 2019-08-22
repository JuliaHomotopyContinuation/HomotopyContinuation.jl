module HomotopyContinuation

    import DynamicPolynomials
    import ElasticArrays
    import FixedPolynomials
    import LinearAlgebra
    import MixedSubdivisions
    import MultivariatePolynomials
    import PrettyTables
    import Printf
    import ProjectiveVectors
    import Random
    import StaticArrays
    import StaticPolynomials
    import TreeViews

    import LinearAlgebra: cond
    import Parameters: @pack!, @unpack
    import DynamicPolynomials: @polyvar, subs, differentiate
    import ProjectiveVectors: PVector
    import StaticArrays: SVector, @SVector
    import Test: @test
    import MixedSubdivisions: mixed_volume

    const FP = FixedPolynomials
    const MP = MultivariatePolynomials
    const SP = StaticPolynomials

    export @polyvar, subs, differentiate
    export mixed_volume
    export cond

    include("progress_meter.jl")
    import .ProgressMeter

    include("utilities.jl")
    include("norms.jl")
    include("linear_algebra.jl")
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
end
