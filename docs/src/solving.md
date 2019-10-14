# Solving general polynomial systems

The [`solve`](@ref) function is the most convenient way to solve general polynomial systems.
For the mathematical background take a look at our [introduction](https://www.juliahomotopycontinuation.org/guides/introduction/) guide.

## Solve
```@docs
solve
Solver
solver
solver_startsolutions
```

## Result

A call to [`solve`](@ref) returns a [`Result`](@ref):
```@docs
Result
seed
```

The nonsingular solutions are obtained as follows.
```@docs
nonsingular
```

The singular solutions are returned by using the following.
```@docs
singular
```

In order to analyse a `Result` we provide the following additional helper functions
```@docs
results
mapresults
solutions
real_solutions
finite
Base.real(::HomotopyContinuation.Results)
at_infinity
failed
multiplicities!(::HomotopyContinuation.Result)
multiplicities(::HomotopyContinuation.Results)
```

If you are interested in the number of solutions of a certain kind we
also provide the following helper functions.
```@docs
nresults
nfinite
nreal
nsingular
nnonsingular
nat_infinity
nfailed
```

Also make sure to check the documentation for [`PathResult`](@ref).

## Estimate the complexity
We provide methods to compute the maximal number of solutions of polynomial systems.
```@docs
bezout_number
mixed_volume
```
