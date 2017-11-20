export AffinePredictorCorrector

"""
    AffinePredictorCorrector

This algorithm uses as an prediction step an explicit Euler method.
Therefore the prediciton step looks like
```math
x_{k+1} = x_k - ΔtJ_H(x, t)^{-1}\\frac{∂H}{∂t}(x, t)
```
and the correction step looks like
```math
x_{k+1} = x_{k+1} - J_H(x, t)^{-1}H(x, t)
```

This algorithm tracks the path in the affine space.
"""
struct AffinePredictorCorrector <: AbstractPathtrackingAlgorithm
end

is_projective(::AffinePredictorCorrector) = false
