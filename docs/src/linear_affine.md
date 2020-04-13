## Coordinates

Linear and affine subspaces can be represented in either extrinsic coordinates ``x`` with ``Ax = b`` or in intrinsic coordinates ``u`` with ``Bu+p``.
To specify which coordinates are given / expected the following singletons can be used:
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
