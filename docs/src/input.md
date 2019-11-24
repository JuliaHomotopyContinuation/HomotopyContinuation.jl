# Input

The [`solve`](@ref) and [`monodromy_solve`](@ref) functions in `HomotopyContinuation.jl`
accept multiple input formats for polynomial systems. These are

* Arrays of polynomials following the [`MultivariatePolynomials`](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface. We export the `@polyvar` macro from the
  [`DynamicPolynomials`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) package to
  create polynomials with the `DynamicPolynomials` implementation.
* Systems constructed with our own symbolic modeling language as implemented in the `ModelKit` module.
  We export the `@var` macro to create variables in this modeling language.
* Systems (resp. homotopies) following the [`AbstractSystem`](@ref) (resp. [`AbstractHomotopy`](@ref))
  interface.


The difference between the `MultivariatePolynomials` and `ModelKit` input
is best shown on an example.
Assume we want to solve the polynomial system

```math
F(x,y) = \begin{bmatrix}
  (2x + 3y + 2)^2 (4x - 2y + 3) \\
  (y - 4x - 5)^3 - 3x^2 + y^2
 \end{bmatrix} = 0
```

Using the `@polyvar` macro from `DynamicPolynomials` we can do
```julia
@polyvar x y
F = [
  (2x + 3y + 2)^2 * (4x - 2y + 3),
  (y - 4x - 5)^3 - 3x^2 + y^2
]
```
```
2-element Array{Polynomial{true,Int64},1}:
 16x³ + 40x²y + 12xy² - 18y³ + 44x² + 68xy + 3y² + 40x + 28y + 12    
 -64x³ + 48x²y - 12xy² + y³ - 243x² + 120xy - 14y² - 300x + 75y - 125
```
We see that our expression got automatically expanded into a monomial basis. Sometimes
this is very useful, but for the fast evaluation of a polynomial system this is not so useful!
Here, `ModelKit` comes into play.
```julia
@var x y
F = [
  (2x + 3y + 2)^2 * (4x - 2y + 3),
  (y - 4x - 5)^3 - 3x^2 + y^2
]
```
```
2-element Array{HomotopyContinuation.ModelKit.Operation,1}:
 (2x + 3y + 2) ^ 2 * ((4x - 2y) + 3)     
 (((y - 4x) - 5) ^ 3 - 3 * x ^ 2) + y ^ 2
```
Compared to the polynomial input we see that it doesn't forget the structure of the input.

For the internal computations both formulations will be converted into efficient straight
line programs. However, from the polynomial formulation we will not be able to recover
the actual formulation of the problem and therefore the generated straight line program
will be less efficient than the one created by `ModelKit`.

However, there are also cases where the polynomial input is preferable. An example is
when you only have an expression of your polynomial system in the monomial basis.
In this case the polynomial input will generate more efficient code since it is more
optimized for this case.

Besides the different macros to generate variables both packages provide a common set
of helpful functions for modeling problems:

* `variables(f, parameters = [])` to obtain a list of all variables.
* `nvariables(f, parameters = [])` to obtain the number of variables.
* `differentiate(f, vars)` to compute the gradient with respect to the given variables
* `subs(f, var => expr)` to substitute variables with expressions
* `monomials(vars, d; homogenous = false)` create all monomials of degree up to `d`
  (resp. exactly degree `d` if `homogenous` = true)

!!! warning "Variable ordering"
    While `MultivariatePolynomials` orders variables in the order of creation, `ModelKit` orders them *alphabetically*. Also in `ModelKit` two variables with the same name are always identical.


## ModelKit

```@docs
ModelKit.@var
ModelKit.@unique_var
ModelKit.System
ModelKit.Homotopy
ModelKit.compile
```
