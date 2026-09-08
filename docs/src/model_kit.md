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

## Supported operations

`Expression`s support the usual arithmetic operations (`+`, `-`, `*` and `/`), powers `x^r` with a constant exponent  `r`, as well as the following elementary functions:

```julia
sin(x)
cos(x)
tan(x)
asin(x)
acos(x)
sinh(x)
cosh(x)
tanh(x)
exp(x)
log(x)
sqrt(x)
```

!!! note "Branch conventions"
    The multivalued complex fuctions `asin`, `acos`, `log` and `sqrt`, as well as non-integer powers, are evaluated using their [principal branches](https://en.wikipedia.org/wiki/Principal_value). For these functions to vary analytically during path tracking, the arguments must stay in the domain of the principal branch. In particular, arguments must not cross a branch cut or pass through a branch point.

## Methods
```@docs
coefficients(f::Expression, vars::AbstractVector{Variable})
coeffs_as_dense_poly
degree(f::Expression, vars::AbstractVector{Variable})
degrees(::AbstractVector{Expression})
differentiate(expr::ModelKit.Basic, vars::AbstractVector{Variable})
Base.conj
dense_poly
evaluate(expr::AbstractArray{<:ModelKit.Basic}, args...)
expand
exponents_coefficients
poly_from_exponents_coefficients
horner
nvariables(::Expression)
monomials(vars::AbstractVector{<:Union{Variable,Expression}}, d::Integer)
subs(ex::ModelKit.Basic, args...)
rand_poly
to_dict
to_number
variables(::Expression)
is_polynomial(::Expression)
get_num_den(::Expression)
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
