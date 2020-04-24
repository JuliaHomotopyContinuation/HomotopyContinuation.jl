# Systems and Homotopies

Systems ([`AbstractSystem`](@ref)) and homotopies ([`AbstractHomotopy`](@ref)) are used for the numerical computations.

## Systems

Systems ([`AbstractSystem`](@ref)) are the basic building blocks of homotopies.

```@docs
AbstractSystem
```

### AffineChartSystem
```@docs
AffineChartSystem
on_affine_chart(F::System, dims)
```

### CompositionSystem
```@docs
CompositionSystem
compose
```

### FixedParameterSystem
```@docs
FixedParameterSystem
fix_parameters
```

### ModelKitSystem
```@docs
ModelKitSystem
```

### RandomizedSystem
```@docs
RandomizedSystem
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
on_affine_chart(F::Homotopy, dims)
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

### ParameterHomotopy
```@docs
ParameterHomotopy
```

### StraightLineHomotopy
```@docs
StraightLineHomotopy
```
