# ModelKit

## Expressions and Variables
```@docs
Expression
Variable
@var
@unique_var
```

## Functions
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

## Index
```@index
Pages   = ["model_kit.md"]
Modules = [ModelKit]
Order   = [:type, :macro, :function]
```
