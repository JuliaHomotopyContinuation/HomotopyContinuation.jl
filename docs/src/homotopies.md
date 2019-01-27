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
ParameterHomotopy
```

We also provide more specialised homotopies, which are mostly used internally currently
but could be useful in conjunction with the [`PathTracker`](@ref) primitive.
```@docs
PatchedHomotopy
```

## Interface for custom homotopies

The great thing is that you are not limited to the homotopies provided by default.
You can define your own homotopy by defining a struct with super type [`AbstractHomotopy`](@ref).
For this the following interface has to be defined.

### Types
```@docs
AbstractHomotopy
AbstractHomotopyCache
NullCache
```

### Mandatory
The following methods are mandatory to implement.
```@docs
cache
evaluate!
jacobian!
dt!
Base.size(::AbstractHomotopy)
```
### Optional
```@docs
evaluate_and_jacobian!
evaluate_and_jacobian
jacobian_and_dt!
evaluate
jacobian
dt
basehomotopy
```
