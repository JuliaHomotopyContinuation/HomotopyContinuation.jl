# Miscellaneous

## Newton's method

```@docs
newton
NewtonResult
is_success(::NewtonResult)
NewtonCache
```

## Norms

```@docs
AbstractNorm
InfNorm
EuclideanNorm
WeightedNorm
distance(u, v, ::AbstractNorm)
norm(u, ::AbstractNorm)
init!(::WeightedNorm, ::AbstractVector)
update!(::WeightedNorm, ::AbstractVector)
```
