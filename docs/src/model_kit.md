# ModelKit

ModelKit is the symbolic input and modeling language of HomotopyContinuation.jl.
It is designed such that you can easily create an efficient formulation of your problem.

## Expressions and Variables
```@docs
Expression
Variable
@var
@unique_var
```

## Methods
```@docs
coefficients
degree
degrees(::AbstractVector{Expression})
differentiate
dense_poly
evaluate
expand
exponents_coefficients
horner
nvariables(::Expression)
monomials
subs
rand_poly
to_dict
to_number
variables(::Expression)
```

## System
```@docs
System
degrees(F::System)
expressions(F::System)
nparameters(::System)
nvariables(::System)
parameters(::System)
variables(::System)
```

## Homotopy
```@docs
Homotopy
expressions(::Homotopy)
nparameters(::Homotopy)
nvariables(::Homotopy)
parameters(::Homotopy)
variables(::Homotopy)
```
