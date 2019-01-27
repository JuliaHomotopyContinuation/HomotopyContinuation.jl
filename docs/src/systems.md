# Polynomial systems

Polynomial systems can be represented in numerous ways in a computer and each
representation has certain tradeoffs. For our purposes the most important thing
is that it is *fast* to evaluate the system. Therefore we automatically convert
an input given by [`DynamicPolynomial`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)s
to another representation more suitable for numerically evaluations.
The default is currently [`FPSystem`](@ref).

## Default systems
We provide the following systems by default.
```@docs
FPSystem
SPSystem
FixedHomotopy
```

## Interface for custom systems

The great thing is that you are not limited to the systems provided by default.
Maybe your polynomial system has a particular structure which you want to use to efficiently
evaluate it. For this you can define your own homotopy by defining a
struct with super type [`AbstractSystem`](@ref).
For this the following interface has to be defined.

### Types
```@docs
AbstractSystem
AbstractSystemCache
NullCache
```

### Mandatory
The following methods are mandatory to implement.
```@docs
cache
evaluate!
evaluate
jacobian!
jacobian
Base.size(::AbstractSystem)
```

### Optional
The following methods are mandatory to implement.
The following are optional to implement but usually you want to define at least
[`cache`](@ref).
```@docs
evaluate_and_jacobian!
evaluate_and_jacobian
```
