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

import MultivariatePolynomials: variables, nvariables, monomials
using MultivariatePolynomials: subs, differentiate

using Base: @propagate_inbounds
using LinearAlgebra: cond
using Parameters: @pack!, @unpack
using DoubleFloats: Double64, ComplexDF64
using DynamicPolynomials: @polyvar
using ProjectiveVectors:
    PVector,
    dims,
    dimension_indices,
    fubini_study,
    affine_chart,
    ×,
    component,
    components,
    combine
using StaticArrays: SVector, @SVector
using Test: @test
using MixedSubdivisions: mixed_volume

const FP = FixedPolynomials
const MP = MultivariatePolynomials
const SP = StaticPolynomials

export @polyvar
export variables, nvariables, subs, differentiate, monomials
export mixed_volume
export cond
export ProjectiveVectors,
    PVector,
    dims,
    dimension_indices,
    affine_chart,
    fubini_study,
    ×,
    component,
    components,
    combine

include("progress_meter.jl")
import .ProgressMeter

include("model_kit.jl")
import .ModelKit
using .ModelKit: @var, @unique_var, System, Homotopy, compile, interpreted
export ModelKit, @var, @unique_var, System, Homotopy, compile, interpreted

include("norms.jl")
include("utilities.jl")
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
include("path_info.jl")

end
