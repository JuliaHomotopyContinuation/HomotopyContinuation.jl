# ModelKit

ModelKit is the symbolic input and modeling language of HomotopyContinuation.jl.
It is designed such that you can easily create an efficient formulation of your problem.

## Expressions and Variables
```@docs
Expression
Variable
@var
@unique_var
variables(prefix::Union{Symbol,String}, indices...)
```

## Methods
```@docs
coefficients(f::Expression, vars::AbstractVector{Variable})
coeffs_as_dense_poly
degree(f::Expression, vars::AbstractVector{Variable})
degrees(::AbstractVector{Expression})
differentiate(expr::ModelKit.Basic, vars::AbstractVector{Variable})
dense_poly
evaluate(expr::AbstractArray{<:ModelKit.Basic}, args...)
expand
exponents_coefficients
horner
nvariables(::Expression)
monomials(vars::AbstractVector{<:Union{Variable,Expression}}, d::Integer)
subs(ex::ModelKit.Basic, args...)
rand_poly
to_dict
to_number
variables(::Expression)
```

## System
```@docs
System
evaluate(F::System, x, p = nothing)
jacobian(F::System)
jacobian(F::System, x, p = nothing)
degrees(F::System)
expressions(F::System)
optimize(::System)
multi_degrees(::System)
nparameters(::System)
nvariables(::System)
parameters(::System)
support_coefficients(::System)
variables(::System)
variable_groups(::System)
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
