export AbstractPredictorCorrectorAlgorithm
abstract type AbstractPredictorCorrectorAlgorithm{Val} <: AbstractHomotopyContinuationAlgorithm end

const APCA{P} = AbstractPredictorCorrectorAlgorithm{P}

include("predictorcorrector/algorithms/interface.jl")
include("predictorcorrector/algorithms/affine.jl")
include("predictorcorrector/algorithms/spherical.jl")
include("predictorcorrector/algorithms/fixed_stepsize.jl")

include("predictorcorrector/trackpath.jl")
include("predictorcorrector/endgames/cauchy.jl")
include("predictorcorrector/result.jl")
include("predictorcorrector/solve.jl")
