# Systems

Systems ([`AbstractSystem`](@ref)) are the basic building blocks of homotopies.

```@docs
AbstractSystem
```

Although they sound similar, [`AbstractSystem`](@ref)s and [`System`](@ref) share different
purposes. [`AbstractSystem`](@ref)s are intented for the fast numerical evaluation
of a fixed system. Whereas a [`System`](@ref) is intended for formulating your problem
symbolically.
A [`System`](@ref) can be converted to two different basic `AbstractSystem`s,
a [`CompiledSystem`](@ref) (fast, but introduce compilation overhead) and an
[`InterpretedSystem`](@ref) (slower, but not compilation overhead).

```@docs
CompiledSystem
InterpretedSystem
MixedSystem
fixed(::System)
```

It is also possible to set the default compile flag.
```@docs
set_default_compile
```

## Interface

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

## AffineChartSystem
```@docs
AffineChartSystem
on_affine_chart(F::System, dims)
```

## CompositionSystem
```@docs
CompositionSystem
compose
```

## FixedParameterSystem
```@docs
FixedParameterSystem
fix_parameters(F::AbstractSystem, p)
```

## RandomizedSystem
```@docs
RandomizedSystem
```

## Real Systems
```@docs
is_real
```
