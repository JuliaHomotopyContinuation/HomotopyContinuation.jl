# Solving Polynomial Systems

At the heart of the package is the [`solve`](@ref) function. It takes
a bunch of different input combinations and returns an [`AffineResult`](@ref) or [`ProjectiveResult`](@ref) depending on the input.

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
solutions
failed
atinfinity
singular
nonsingular
finite
seed
```

If you are interested in the number of solutions of a certain kind we
also provide the following helper functions.
```@docs
nresults
nfinite
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
issuccess
isfailed
isaffine
isprojective
isatinfinity
issingular
isnonsingular
```
