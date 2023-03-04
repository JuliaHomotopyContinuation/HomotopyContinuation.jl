# Positive Dimensional Solution Sets

The basic data structure to work with positive dimensional solution sets ``V(F)`` of a polynomial system ``F`` is a [`WitnessSet`](@ref) ``W``.
The general idea is to intersect ``V(F)`` with an (affine) linear space ``L`` such that
the intersection ``V(F) ∩ L`` consists of only finitely many points (witnesses).

```@docs
WitnessSet
```

To compute a [`WitnessSet`](@ref) call [`witness_set`](@ref). 

```@docs
witness_set
```

To compute witness sets for all possible dimensions use [`regeneration`](@ref) or [`witness_sets`](@ref).

```@docs
regeneration
witness_sets
```

Witness sets can be decomposed into [irreducible components](https://en.wikipedia.org/wiki/Irreducible_component) by using the [`decompose`](@ref) function.

```@docs
decompose
```

## Numerical Irreducible Decomposition

A numerical irreducible decomposition for ``F`` is a collection of witness sets ``W₁, ..., Wᵣ`` such that the ``Wᵢ`` are all witness sets for different [irreducible components](https://en.wikipedia.org/wiki/Irreducible_component) and ``V(F)`` is their union. 

```@docs
numerical_irreducible_decomposition
nid
```

The output of [`nid`](@ref) is a [`NumericalIrreducibleDecomposition`](@ref).

```@docs
NumericalIrreducibleDecomposition
```

```@docs
n_components
witness_sets
degrees
```


## Computing Information about Witness Sets

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
