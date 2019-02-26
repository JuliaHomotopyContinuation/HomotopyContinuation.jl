# Reference

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

## Utilities

```@docs
bezout_number
ishomogenous
uniquevar
homogenize
```

## AffinePatches
Affine patches are there to augment projective system such that they can be considered
as (locally) affine system. By default the following patches are defined

```@docs
OrthogonalPatch
EmbeddingPatch
RandomPatch
FixedPatch
```
