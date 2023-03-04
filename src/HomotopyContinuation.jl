module HomotopyContinuation

export ModelKit

using Arblib: Arblib
import DelimitedFiles
using DynamicPolynomials: @polyvar
import ElasticArrays: ElasticArrays, ElasticArray
import FiniteDiff
import LinearAlgebra
import LoopVectorization
using LRUCache: LRU
import MixedSubdivisions: MixedSubdivisions, MixedCell, mixed_volume
import MultivariatePolynomials
const MP = MultivariatePolynomials
using Parameters: @unpack
import ProgressMeter
const PM = ProgressMeter
import Random
import Printf
import PrettyTables
using Reexport: @reexport
import SemialgebraicSets
import StructArrays
import TreeViews
import Base.Iterators: product, reverse
import Base: append!

const LA = LinearAlgebra

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
include("voronoi_tree.jl")
include("unique_points.jl")
include("linear_algebra.jl")
include("linear_subspaces.jl")
include("systems.jl")
include("homotopies.jl")
include("predictor.jl")
include("newton_corrector.jl")
include("newton.jl")
include("tracker.jl")
include("valuation.jl")
include("path_result.jl")
include("endgame_tracker.jl")
include("path_info.jl")
include("total_degree.jl")
include("binomial_system.jl")
include("polyhedral.jl")
include("overdetermined.jl")
include("result.jl")
include("solve.jl")
include("monodromy.jl")
include("witness_set.jl")
include("certification.jl")
include("semialgebraic_sets.jl")
include("numerical_irreducible_decomposition.jl")

function __init__()
    # Overwrite default IJulia behaviour. See
    # https://github.com/timholy/ProgressMeter.jl/issues/162
    if isdefined(ProgressMeter, :ijulia_behavior)
        ProgressMeter.ijulia_behavior(:clear)
    end
end

end #module
