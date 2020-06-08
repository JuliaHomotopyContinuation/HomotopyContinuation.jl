module HomotopyContinuation

export ModelKit

using DynamicPolynomials: @polyvar
import ElasticArrays: ElasticArray
import FiniteDiff
import LinearAlgebra
import LoopVectorization
import MixedSubdivisions
import MultivariatePolynomials
const MP = MultivariatePolynomials
using Parameters: @unpack
import ProgressMeter
import Random
import Printf
import PrettyTables
using Reexport: @reexport
import StructArrays
import TreeViews

const LA = LinearAlgebra

include("ModelKit.jl")
export @polyvar

include("DoubleDouble.jl")
using .DoubleDouble

using ProjectiveVectors: PVector, dims, dimension_indices

include("utils.jl")
include("norm.jl")
include("voronoi_tree.jl")
include("unique_points.jl")
include("linear_algebra.jl")
include("linear.jl")
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
include("solve.jl")
include("result.jl")
include("monodromy.jl")

end
