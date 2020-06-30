# Homotopies

Homotopies ([`AbstractHomotopy`](@ref)) are at the heart of homotopy continuation.
A homotopy is a parameterized family ``H(x,t)`` of polynomial systems.
By convention, homotopies are tracked from ``t=1`` to ``t=0``, i.e., ``H(x,1)`` is considered
the *start system* and ``H(x,0)`` is the *target system*.

```@docs
AbstractHomotopy
```
As for [`AbstractSystem`](@ref)s, [`AbstractHomotopy`](@ref) and [`Homotopy`](@ref) share different purposes.
A [`Homotopy`](@ref) is intended for formulating your problem
symbolically.
A [`Homotopy`](@ref) can be converted to two different basic `AbstractHomotopy`s,
a [`CompiledHomotopy`](@ref) (fast, but introduce compilation overhead) and an
[`InterpretedHomotopy`](@ref) (slower, but not compilation overhead).

```@docs
CompiledHomotopy
InterpretedHomotopy
fixed(::Homotopy)
```

## AffineChartHomotopy
```@docs
AffineChartHomotopy
on_affine_chart(F::Homotopy, dims)
```

## CoefficientHomotopy
```@docs
CoefficientHomotopy
```

## FixedParameterHomotopy
```@docs
FixedParameterHomotopy
fix_parameters(H::AbstractHomotopy, p)
```

## IntrinsicSubspaceHomotopy
```@docs
IntrinsicSubspaceHomotopy
set_subspaces!
```

## ParameterHomotopy
```@docs
ParameterHomotopy
```

## StraightLineHomotopy
```@docs
StraightLineHomotopy
```
