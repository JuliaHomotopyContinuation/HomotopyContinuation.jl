# Solving Polynomial Systems


## solve
```@docs
solve
```


## Result

A call to [`solve`](@ref) returns a [`Result`](@ref):
```@docs
Result
seed
```

In order to analyse a `Result` we provide the following helper functions
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

Also make sure to check the documentation for [`PathResult`](@ref).
