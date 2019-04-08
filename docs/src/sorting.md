# Sorting arrays of solutions

We provide functions for sorting analyzing arrays of vectors.

## Computing unique points in an array of vectors
```@docs
UniquePoints
```
We provide several helper functions for `UniquePoints`.
```@docs
points
iscontained
add!
simple_add!
empty!
```

## Computing points in an array of vectors which appear multiple times
If instead of unique points, the user wants to have the information which points in an array of points appear with multiplicity, they should use the next function.
```@docs
multiplicities
```
The `multiplicities` functions may also be applied to [`Result`](@ref); see here:
[`multiplicities(::HomotopyContinuation.Results)`](@ref).
