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

We also provide more specialised homotopies, which are mostly used internally currently
but could be useful in conjunction with the [`PathTracking.PathTracker`](@ref) primitive.
```@docs
Homotopies.PatchedHomotopy
Homotopies.PatchSwitcherHomotopy
```

## Interface for custom homotopies

The great thing is that you are not limited to the homotopies provided by default.
You can define your own homotopy by defining a struct with super type [`Homotopies.AbstractHomotopy`](@ref).
For this the following interface has to be defined.

### Types
```@docs
Homotopies.AbstractHomotopy
Homotopies.AbstractHomotopyCache
Homotopies.NullCache
```

### Mandatory
The following methods are mandatory to implement.
```@docs
Homotopies.cache
Homotopies.evaluate!
Homotopies.jacobian!
Homotopies.dt!
Base.size(::Homotopies.AbstractHomotopy)
```
### Optional
```@docs
Homotopies.evaluate_and_jacobian!
Homotopies.evaluate_and_jacobian
Homotopies.jacobian_and_dt!
Homotopies.evaluate
Homotopies.jacobian
Homotopies.dt
Homotopies.precondition!
Homotopies.update!
```
