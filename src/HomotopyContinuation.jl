module HomotopyContinuation

    import DoubleFloats
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

    using Base: @propagate_inbounds
    using LinearAlgebra: cond
    using Parameters: @pack!, @unpack
    using DoubleFloats: Double64, ComplexDF64
    using DynamicPolynomials: @polyvar, subs, differentiate
    using ProjectiveVectors: PVector
    using StaticArrays: SVector, @SVector
    using Test: @test
    using MixedSubdivisions: mixed_volume

    const FP = FixedPolynomials
    const MP = MultivariatePolynomials
    const SP = StaticPolynomials

    export @polyvar, subs, differentiate
    export mixed_volume
    export cond
    export PVector

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

    include("predictors.jl")
    include("newton_corrector.jl")
    include("core_tracker.jl")

    include("valuation.jl")
    include("cauchy_endgame.jl")
    include("path_tracker.jl")
    include("result.jl")
    include("polyhedral.jl")
    include("overdetermined.jl")
    include("solver.jl")
    include("monodromy.jl")
end
