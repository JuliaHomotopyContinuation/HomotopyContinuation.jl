module HomotopyContinuation

export ModelKit

using Arblib: Arblib
using DynamicPolynomials: @polyvar
import ElasticArrays: ElasticArrays, ElasticArray
import LinearAlgebra
import LoopVectorization
import MixedSubdivisions: MixedSubdivisions, MixedCell, mixed_volume
import MultivariatePolynomials
const MP = MultivariatePolynomials
using Parameters: @unpack
import ProgressMeter
const PM = ProgressMeter
import Random
import Printf
using Reexport: @reexport
import StructArrays
import Base.Iterators: product, reverse
import Base: push!

const LA = LinearAlgebra

# To ensure that ModelKit and HomotopyContinuation define methods of the same `is_real` function:
function is_real end

include("DoubleDouble.jl")
using .DoubleDouble

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
include("systems.jl")
include("homotopies.jl")
include("predictor.jl")
include("newton_corrector.jl")
include("newton.jl")
include("tracker.jl")
include("valuation.jl")
include("path_result.jl")
include("endgame_tracker.jl")
include("total_degree.jl")
include("binomial_system.jl")
include("polyhedral.jl")
include("result.jl")
include("solve.jl")

function __init__()
    # Overwrite default IJulia behaviour. See
    # https://github.com/timholy/ProgressMeter.jl/issues/162
    if isdefined(ProgressMeter, :ijulia_behavior)
        ProgressMeter.ijulia_behavior(:clear)
    end
end

end #module
