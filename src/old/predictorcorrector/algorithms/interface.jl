export predict!, correct!

"""
    precondition!(x, H, algorithm)

Precondition the startsolution `x`.
"""
function precondition! end

"""
```
predict!(
    u::Vector{T},
    A::Matrix{T},
    b::Vector{T},
    H::AbstractHomotopy{T},
    J_H!::Function,
    Hdt!::Function,
    x::Vector{T},
    t::Number,
    Δt::Number,
    alg::AbstractPredictorCorrectorAlgorithm)
```

Make a prediction step for the algorithm `alg` and the homotopy `H` with the jacobian
`J_H`, the derivative w.r.t t `∂H∂t` at `x` to the time `t` with a stepwidth `Δt`.
"""
function predict! end

"""
```
correct!(
    u::Vector{T},
    A::Matrix{T},
    b::Vector{T},
    H::AbstractHomotopy,
    J_H!::Function,
    x::Vector{T},
    t::Number,
    tol::Float64,
    max_iterations::Int,
    alg::AbstractPredictorCorrectorAlgorithm)
```
Make a correction step for the algorithm `alg` and the homotopy `H` with the jacobian
`J_H` at `x` to the time `t`. Stores the result in `u` and returns a boolean indicating
whether ``|u-x|<tol`` withing `max_iterations` iterations.
"""
function correct! end
