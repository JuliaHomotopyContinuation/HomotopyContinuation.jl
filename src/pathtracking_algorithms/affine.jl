export AffinePredictorCorrector

"""
    AffinePredictorCorrector

A prediction-correction based algorithm. As prediction it used single step of the explicit
Euler method.
"""
struct AffinePredictorCorrector <: AbstractPathtrackingAlgorithm
end

is_projective(::AffinePredictorCorrector) = false
