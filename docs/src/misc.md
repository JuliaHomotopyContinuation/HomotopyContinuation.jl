# Miscellaneous

## Storing solutions and parameters

```@docs
write_solutions
read_solutions
write_parameters
read_parameters
```

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

## Unique points, group actions and multiplicities

```@docs
UniquePoints
search_in_radius(::UniquePoints, v, tol::Real)
add!(UP::UniquePoints{T,Id,M,GA}, v, id::Id, tol::Real) where {T,Id,M,GA}
multiplicities
unique_points
```

## Debugging

```@docs
path_info
```
