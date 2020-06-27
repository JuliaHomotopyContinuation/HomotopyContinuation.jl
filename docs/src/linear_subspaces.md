# Linear Subspaces

We provide built-in data structures to work with (affine) linear subspaces ``L``
where ``L`` can be represented in either *extrinsic coordinates* ``x`` with
``L = \{ x | Ax = b \}`` or in *intrinsic coordinates* ``u`` with ``L = \{ Bu+p | u \}``.

## Coordinates

To specify which coordinates are given / expected the following can be used:
```@docs   
Coordinates
Intrinsic
Extrinsic
```

## Constructors

```@docs
LinearSubspace
AffineExtrinsic
AffineIntrinsic
```

### Functions

```@docs
ambient_dim
codim
coord_change
dim
intrinsic
is_linear
extrinsic
geodesic
geodesic_distance
rand_subspace
translate
```
