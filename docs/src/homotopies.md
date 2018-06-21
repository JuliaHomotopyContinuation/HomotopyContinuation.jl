# Homotopies

A homotopy is a function
```math
H: \mathbb{C}^N × \mathbb{C} → \mathbb{C}^n, (x,t) ↦ H(x,t)
```
where ``H(⋅, t)`` is a polynomial system for all ``t∈\mathbb{C}``.

## Default homotopies
The following homotopies are available by default
```@docs
StraightLineHomotopy
FixedPointHomotopy
```

## Homotopy interface

The great thing is that you are not limited to the homotopies provided by default.
You can define your own homotopy by defining a struct with super type [`Homotopies.AbstractHomotopy`](@ref).
For this the following interface has to be defined.

### Abstract types
```@docs
Homotopies.AbstractHomotopy
Homotopies.AbstractHomotopyCache
```

### Mandatory
The following methods are mandatory to implement.
```@docs
Homotopies.evaluate!
Homotopies.jacobian!
Homotopies.dt!
Base.size(::Homotopies.AbstractHomotopy)
```
### Optional
The following are optional to implement but usually you want to define at least
[`Homotopies.cache`](@ref).
```@docs
Homotopies.cache
Homotopies.evaluate_and_jacobian!
Homotopies.evaluate_and_jacobian
Homotopies.jacobian_and_dt!
Homotopies.evaluate
Homotopies.jacobian
Homotopies.dt
Homotopies.precondition!
Homotopies.update!
```
