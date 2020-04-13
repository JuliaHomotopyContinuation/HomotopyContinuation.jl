# Linear and Affine Subspaces

We provide built-in data structures to work with affine and linear subspaces ``L``.
``L`` can be represented in either *extrinsic coordinates* ``x`` with
``L =\{ x | Ax = b \}`` or in *intrinsic coordinates* ``u`` with ``L=\{Bu+p | u\}``.

## Coordinates

To specify which coordinates are given / expected the following can be used:
```@docs   
Coordinates
Intrinsic
Extrinsic
```

## Affine Subspace

```@docs
AffineSubspace
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
extrinsic
geodesic
geodesic_distance
rand_affine_subspace
```
