# HomotopyContinuation.jl [WIP]
| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] |


`HomotopyContinuation.jl` is a package for solving square polynomial systems via homotopy continuation.

The aim of this project is twofold: establishing a fast numerical polynomial solver in `Julia` and at the same time providing a highly customizable algorithmic environment for researchers for designing and performing individual experiments.

You can simply install this package via the Julia package manager
```julia
Pkg.add("JuliaHomotopyContinuation");
```

## A first example
HomotopyContinuation.jl aims at having easy-to-understand top-level commands. For instance, suppose we wanted to solve the following system

```math
f= \begin{bmatrix} x^2+y\\ y^2-1\end{bmatrix}.  
```

First, we have to define ``f`` in Julia. For this purpose
HomotopyContinuation.jl provides an interface to [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) for human-readable input and output.

```julia
import DynamicPolynomials: @polyvar # @polyvar is a function for initializing variables.

@polyvar x y # initialize the variables x y
f = [x^2+y, y^2-1]
```

To solve  ``f=0`` we execute the following command.

```julia
using HomotopyContinuation # load the module HomotopyContinuation

solve(f) # solves for f=0
```

(see [here](@ref solveroptions) for a list of options that can be passed to `solve`).

The last command will return a type `HomotopyContinuation.Result{Complex{Float64}}` of length 4 (one entry for each solution):

```julia-repl
julia> ans

julia> HomotopyContinuation.Result{Complex{Float64}}
# paths → 4
# successfull paths → 4
# solutions at infinity → 0
# singular solutions → 0
# real solutions → 2
HomotopyContinuation.PathResult{Complex{Float64}}[4]
```

Let us see what is the information that we get. Four paths were attempted to be solved, four of which were completed successfully. Since we tried to solve an affine system, the algorithm checks whether there are solutions at infinity: in this case there are none. None of the solutions is singular and two of them are real. To access the first solution in the array we write

```julia-repl
julia> ans[1]

julia> HomotopyContinuation.PathResult{Complex{Float64}}
returncode → :success
solution → Complex{Float64}[2]
singular → false
residual → 1.02e-15…
newton_residual → 8.95e-16…
log10_condition_number → 0.133…
windingnumber → 1
angle_to_infinity → 0.615…
real_solution → true
startvalue → Complex{Float64}[2]
iterations → 17
endgame_iterations → 5
npredictions → 2
predictions → Vector{Complex{Float64}}[2]
```

The returncode tells us that the pathtracking was successfull. What do the entries of that table tell us? Let us consider the most relevant (for a complete list of explanations consider [this](@ref result) section).

- `solution`: the zero that is computed (here it is ``[-1,-1]``).
- `singular`: boolean that shows whether the zero is singular.
- `residual`: the computed value of ``|f([-1,-1])|``.
- `angle_to_infinity`: the algorithms homogenizes the system ``f`` and then computes all solutions in projective space. The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
-  `real_solution`: boolean that shows whether the zero is real.

Suppose we were only interested in the real solutions. The command to extract them is

```julia
solutions(solve(f), only_real=true)
```
(a detailed explanation of the `solutions` function is [here](@ref solutions)). Indeed, we have

```julia-repl
julia> [ans[i].solution for i=1:2]
julia> Vector{Complex{Float64}}[2]
Complex{Float64}[2]
1.00… - 2.66e-15…im
-1.00… + 1.33e-15…im
Complex{Float64}[2]
-1.00… + 2.72e-15…im
-1.00… + 1.44e-15…im
```
which are the two real zeros of `f`. By assigning the boolean values in the [`solutions`](@ref solutions) function we can filter the solutions given by `solve(f)` according to our needs.

We solve some more elaborate systems in the [example section](@ref examples).

`JuliaHomotopyContinuation` also supports input of type `BigFloat`.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://JuliaHomotopyContinuation.github.io/HomotopyContinuation.jl/stable
[docs-latest-url]: https://JuliaHomotopyContinuation.github.io/HomotopyContinuation.jl/latest

[build-img]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl
[codecov-img]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl
