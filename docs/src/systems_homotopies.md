# Systems and Homotopies

Systems ([`AbstractSystem`](@ref)) and homotopies ([`AbstractHomotopy`](@ref)) are used for the numerical computations.

## Systems

Systems ([`AbstractSystem`](@ref)) are the basic building blocks of homotopies.

```@docs
AbstractSystem
```
### ModelKitSystem
```@docs
ModelKitSystem
```
### RandomizedSystem
```@docs
RandomizedSystem
excess_solution_check!
excess_solution_check
```


## Homotopies
Homotopies ([`AbstractHomotopy`](@ref)) are at the heart of homotopy continuation.
A homotopy is a parameterized family ``H(x,t)`` of polynomial systems.
By convention, homotopies are tracked from ``t=1`` to ``t=0``, i.e., ``H(x,1)`` is considered
the *start system* and ``H(x,0)`` is the *target system*.

```@docs
AbstractHomotopy
```

### AffineChartHomotopy
```@docs
AffineChartHomotopy
on_affine_chart
```

### AffineSubspaceHomotopy
```@docs
AffineSubspaceHomotopy
set_subspaces!
```

### ModelKitHomotopy
```@docs
ModelKitHomotopy
```

### StraightLineHomotopy
```@docs
ParameterHomotopy
```

### StraightLineHomotopy
```@docs
StraightLineHomotopy
```
