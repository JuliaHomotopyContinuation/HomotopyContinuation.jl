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
WeightedNorm
distance(u, v, ::AbstractNorm)
norm(u, ::AbstractNorm)
init!(::WeightedNorm, ::AbstractVector)
update!(::WeightedNorm, ::AbstractVector)
```

## Debugging

```@docs
path_info
```
