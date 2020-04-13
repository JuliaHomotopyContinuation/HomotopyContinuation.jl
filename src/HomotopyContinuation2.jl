module HomotopyContinuation2

export ModelKit

import ElasticArrays: ElasticArray
import LinearAlgebra
import LoopVectorization
import MixedSubdivisions
using Parameters: @unpack
import Printf
import PrettyTables
using Reexport: @reexport
import StructArrays

const LA = LinearAlgebra

include("ModelKit.jl")
include("DoubleDouble.jl")
using .DoubleDouble

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

include("utils.jl")
include("norm.jl")
include("linear_algebra.jl")
include("linear.jl")
include("systems.jl")
include("homotopies.jl")
include("total_degree.jl")
include("predictor.jl")
include("newton_corrector.jl")
include("newton.jl")
include("tracker.jl")
include("valuation.jl")
include("path_tracker.jl")
include("path_info.jl")
include("binomial_system.jl")
include("polyhedral.jl")


end
