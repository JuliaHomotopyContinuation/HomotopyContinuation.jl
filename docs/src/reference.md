# Reference

## AffinePatches
Affine patches are there to augment projective system such that they can be considered
as (locally) affine system.

```@docs
AbstractAffinePatch
OrthogonalPatch
EmbeddingPatch
RandomPatch
FixedPatch
```

## Input
We support any polynomials which follow the [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl)
interface. By default we export the routines `@polyvar`, `PolyVar`, `differentiate`
and `variables`
from the [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)
implementation.
With these you can simply create variables
```julia
# Create variables x, y, z
@polyvar x y z
f = x^2+y^2+z^2

# You can also create an array of variables
@polyvar x[1:3] # This creates x1, x2, x3 accessed by x[1], x[2], x[3]
f = dot(x, x) # = x[1]^2+x[2]^2+x[3]^2

# Also you can create matrices of variables
# This creates x1_1, x1_2, x2_1, x2_2 accessed by
# x[1,1], x[1,2], x[2,1], x[2,2]
@polyvar x[1:2, 1:2]
```

We also provide methods construct compositions of polynomial systems:
```@docs
compose
```

## Distances and norms

We provide functions for computing norms and distance.
These are subtypes of `AbstractNorm`:

```@docs
AbstractNorm
```

They implement
```@docs
distance(u, v, norm::AbstractNorm)
LinearAlgebra.norm(x, ::AbstractNorm)
```

The following norms are implemented:
```@docs
EuclideanNorm
InfNorm
WeightedNorm
weights(WN::WeightedNorm)
init!(w::WeightedNorm, x::AbstractVector)
update!(w::WeightedNorm, x::AbstractVector)
```

### Deprecated

```@docs
euclidean_distance
euclidean_norm
```

## Path informations

```@docs
path_info
```

## Polynomial Utilities

```@docs
homogenize
is_homogeneous
uniquevar
linear_system
```

## Predictors and Correctors

We use a predictor-corrector scheme to track paths. While we have a fixed implementation of Newton's method as a corrector there multiple predictors available:

```@docs
Euler
Heun
Ralston
RK3
RK4
Pade21
NullPredictor
```

## Projective vectors

To represent solutions in projective spaces we use the
[ProjectiveVectors.jl](https://github.com/JuliaHomotopyContinuation/ProjectiveVectors.jl)
package. The following methods are re-exported from the package.

```@docs
PVector
dims
dimension_indices
affine_chart
×
component
components
combine
```

## Unique points and multiplicities

We provide functions for sorting analyzing arrays of vectors.

### Unique points

We provide the [`unique_points`](@ref) methods to sort through solutions.
```@docs
unique_points
```

The [`unique_points`](@ref) method is powered by the [`UniquePoints`](@ref) data structure.

```@docs
UniquePoints
```

We provide several helper functions for `UniquePoints`.
```@docs
points
is_contained
add!
simple_add!
empty!
```

### Multiplicities

If instead of unique points, the user wants to have the information which points in an array of points appear with multiplicity, they should use the next function.
```@docs
multiplicities
```
The `multiplicities` functions may also be applied to [`Result`](@ref); see here:
[`multiplicities(::HomotopyContinuation.Results)`](@ref).
