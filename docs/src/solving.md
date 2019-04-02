# Solving Polynomial Systems

At the heart of the package is the [`solve`](@ref) function. It takes
a bunch of different input combinations and returns an [`AffineResult`](@ref) or [`ProjectiveResult`](@ref) depending on the input.

The `solve` function works in 3 stages.

1) It takes the input and constructs a homotopy ``H(x,t)`` such that ``H(x,1)=G(x)`` and ``H(x,0)=F(x)`` as well as start solutions ``\mathcal{X}`` where for all ``x_1 ∈ \mathcal{X}`` we have ``H(x_1, 1) ≈ 0``. This step highly depends on the concrete input you provide.

2) Then all start solutions ``x(1) = x_1 ∈ \mathcal{X} `` will be tracked to solutions ``x({t_e})`` such that ``H(x({t_e}), t_e) ≈ 0`` using a predictor-corrector scheme where ``t_e`` is a value between ``0`` and ``1`` (by default ``0.1``).

3) From these intermediate solutions ``x({t_e})`` we start the so called *endgame*. This is an algorithm to predict the value ``x(0)`` for each intermediate solution.

The reason for step 3 is that the final solutions can be singular which provides significant challenges for our gradient based predictor-corrector methods. For more background also check the [FAQ](http://localhost:1313/faq/).

## The *solve* function
```@docs
solve
```

## The result of *solve*

Depending on the input `solve` returns one of the following types
```@docs
AffineResult
ProjectiveResult
```
A `Result` is a wrapper around the results of each single path ([`PathResult`](@ref)) and it contains some additional informations like
the used random seed for the computation.

In order to analyze a `Result` we provide the following helper functions
```@docs
results
mapresults
solutions
realsolutions
uniquesolutions
finite
Base.real(::HomotopyContinuation.Results)
atinfinity
singular
nonsingular
failed
multiplicities(::HomotopyContinuation.Results)
seed
```

If you are interested in the number of solutions of a certain kind we
also provide the following helper functions.
```@docs
nresults
nfinite
nreal
nsingular
nnonsingular
natinfinity
nfailed
```

### PathResult
For each path we return a [`PathResult`](@ref) containing the detailed information about
the single path.
```@docs
PathResult
```

The following helper functions are provided
```@docs
solution
residual
startsolution
Base.isreal(::PathResult)
LinearAlgebra.issuccess(::PathResult)
isfailed
isaffine
isprojective
isatinfinity
issingular
isnonsingular
```
