export SphericalPredictorCorrector

"""
    SphericalPredictorCorrector


A prediction-correction based algorithm. As prediction it uses explicit euler.
Jacobian is augmented by ``x{^H}`` to solve a square system during correction.
**MORE DOCUMENTATION NEEDED**
"""
struct SphericalPredictorCorrector <: AbstractPathtrackingAlgorithm
end

is_projective(::SphericalPredictorCorrector) = true
