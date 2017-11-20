export SphericalPredictorCorrector

"""
    SphericalPredictorCorrector

This algorithm uses as an prediction step an explicit Euler method.
For the prediction and correction step the Jacobian is augmented by Hermitian transposed
``x^*`` to get a square system.
Therefore the prediciton step looks like
```math
x_{k+1} = x_k - Δt\\begin{bmatrix}
    J_H(x, t) \\\\
    x^*
\\end{bmatrix}^{-1}
\\frac{∂H}{∂t}(x, t)
```
and the correction step looks like
```math
x_{k+1} = x_{k+1} - \\begin{bmatrix}
    J_H(x, t) \\\\
    x^*
\\end{bmatrix}^{-1}
H(x, t)
```

After each prediciton and correction the algorithm normalizes ``x`` again, i.e. ``x`` is
then a point on the sphere again.

This algorithm tracks the path in the projective space.
"""
struct SphericalPredictorCorrector <: AbstractPathtrackingAlgorithm
end

is_projective(::SphericalPredictorCorrector) = true
