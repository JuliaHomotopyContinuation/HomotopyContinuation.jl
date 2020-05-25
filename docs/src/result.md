# Results

A call to [`solve`](@ref) returns a [`Result`](@ref):
```@docs
Result
seed(::Result)
path_results(::Result)
```

## Filtering results and solutions
```@docs
results
solutions
Base.real(::Result)
real_solutions
nonsingular
singular
at_infinity
failed
```

## Counting
```@docs
nresults
nsolutions
nreal
nnonsingular
nsingular
nat_infinity
nexcess_solutions
nfailed
```

## PathResult

```@docs
PathResult
solution(::PathResult)
is_success(::PathResult)
is_at_infinity(::PathResult)
is_excess_solution(::PathResult)
is_failed(::PathResult)
is_finite(::PathResult)
is_singular(::PathResult)
is_nonsingular(::PathResult)
is_real(::PathResult)
accuracy(::PathResult)
residual(::PathResult)
steps(::PathResult)
accepted_steps(::PathResult)
rejected_steps(::PathResult)
winding_number(::PathResult)
path_number(::PathResult)
start_solution(::PathResult)
multiplicity(::PathResult)
last_path_point(::PathResult)
valuation(::PathResult)
```
