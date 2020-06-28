# Witness Sets

A [`WitnessSet`](@ref) ``W`` is the basic data structure to work with positive dimensional solution sets ``V(F)`` of a polynomial system ``F``.
The general idea is to intersect ``V(F)`` with an (affine) linear space ``L`` such that
the intersection ``V(F) ∩ L`` consists of only finitely many points (witnesses).
Over the complex numbers the number of points is independent of the linear space ``L`` and called the *degree* of ``V(F)``.

```@docs
WitnessSet
```

To compute a [`WitnessSet`](@ref) call [`witness_set`](@ref).

```@docs
witness_set
```

To obtain information about a [`WitnessSet`](@ref) the following functions are provided.
```@docs
solutions(::WitnessSet)
results(W::WitnessSet{<:Any,<:Any,PathResult})
system
linear_subspace
dim(::WitnessSet)
codim(::WitnessSet)
ModelKit.degree(::WitnessSet)
```

To test for completeness of a [`WitnessSet`](@ref) you can perform a [`trace_test`](@ref)
```@docs
trace_test
```
