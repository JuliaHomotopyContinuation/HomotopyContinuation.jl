module HomotopyContinuation

export ModelKit


import DelimitedFiles
import FiniteDiff
import IntervalTrees
import LinearAlgebra
import LoopVectorization
import MultivariatePolynomials
import PrettyTables
import Printf
import ProgressMeter
import Random
import SemialgebraicSets
import StaticArrays
import StructArrays
import TreeViews

using Arblib: Arblib
using Base: RefValue, push!
using Base.Iterators: product, reverse
using DynamicPolynomials: @polyvar
using ElasticArrays: ElasticArrays, ElasticArray
using LRUCache: LRU
using MixedSubdivisions: MixedSubdivisions, MixedCell, mixed_volume
using Parameters: @unpack
using Reexport: @reexport

const LA = LinearAlgebra
const MP = MultivariatePolynomials
const PM = ProgressMeter

include("DoubleDouble.jl")
using .DoubleDouble

include("interval_arithmetic.jl")
using .IntervalArithmetic

include("ModelKit.jl")
export @polyvar

using ProjectiveVectors: PVector, dims, dimension_indices

const COMPILE_DEFAULT = Ref(:mixed)

export set_default_compile

"""
    set_default_compile(mode::Symbol)

Set the default value for the `compile` flag in [`solve`](@ref) and other functions.
Possible values are `:mixed` (default), `:all` and `:none`.
"""
function set_default_compile(mode::Symbol)
    mode ∈ [:mixed, :all, :none] ||
        error("Invalid value `:$mode`, valid values are `:mixed`, `:all`, `:none`.")
    COMPILE_DEFAULT[] = mode
end

include("utils.jl")
include("norm.jl")
include("linear_algebra.jl")
include("linear_solve_workspace.jl")
include("linear_subspaces.jl")
include("systems.jl")
include("homotopies.jl")
include("voronoi_tree.jl")
include("unique_points.jl")
# include("newton.jl")
include("path_result.jl")
include("tracking/tracker.jl")
include("endgame/endgame.jl")
include("total_degree.jl")



include("binomial_system.jl")
include("polyhedral.jl")
# include("overdetermined.jl")
# include("result.jl")
# include("solve.jl")
# include("monodromy.jl")
# include("witness_set.jl")
# include("certification.jl")
# include("semialgebraic_sets.jl")
# include("numerical_irreducible_decomposition.jl")

function __init__()
    # Overwrite default IJulia behaviour. See
    # https://github.com/timholy/ProgressMeter.jl/issues/162
    if isdefined(ProgressMeter, :ijulia_behavior)
        ProgressMeter.ijulia_behavior(:clear)
    end
end

end #module
