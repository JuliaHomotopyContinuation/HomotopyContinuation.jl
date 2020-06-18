# Systems and Homotopies

Systems ([`AbstractSystem`](@ref)) and homotopies ([`AbstractHomotopy`](@ref)) are used for the numerical computations.

## Systems

Systems ([`AbstractSystem`](@ref)) are the basic building blocks of homotopies.

```@docs
AbstractSystem
```

### Interface

An [`AbstractSystem`](@ref) needs to implement the following methods:

```
Base.size(F::AbstractSystem)
ModelKit.variables(F::AbstractSystem)::Vector{Variable}
ModelKit.parameters(F::AbstractSystem) = Variable[]
ModelKit.variable_groups(::AbstractSystem)::Union{Nothing,Vector{Vector{Variable}}} = nothing
 # this has to work with x::Vector{Variable}
(F::AbstractSystem)(x, p = nothing)
 # this has to work with x::Vector{ComplexF64} and x::Vector{ComplexDF64}
evaluate!(u, F::AbstractSystem, x, p = nothing)
# this has only to work with x::Vector{ComplexF64}
evaluate_and_jacobian!(u, U, F::AbstractSystem, x, p = nothing)
```

If the system should be used in context of a parameter homotopy it is also necessary to
implement

```
taylor!(u, ::Val{1}, F::AbstractSystem, x, p::TaylorVector{2})
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
