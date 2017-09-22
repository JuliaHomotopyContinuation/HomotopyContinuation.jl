# HomotopyContinuation.jl [WIP]
[![Build Status](https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl.svg?branch=totaldegree)](https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl)
[![codecov](https://codecov.io/gh/JuliaHomotopyContinuation/HomotopyContinuation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaHomotopyContinuation/HomotopyContinuation.jl)

This is a package for solving systems of polynomials via homotopy continuation.

Getting started
-----------
Since this package is pre-release and also relies on couple of unreleased packages. To satisfy all dependencies you have to install it via
```sh
Pkg.clone("https://github.com/JuliaHomotopyContinuation/Homotopy.jl");
Pkg.clone("https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git")
```

And here is a simple example:
```julia
using HomotopyContinuation
using DynamicPolynomials

@polyvar x
f = (x - 2.0) * (x - (2.5+ 4.0im))
g = (x - 4.3im) * (x + (2.1 - 4im))

H = StraightLineHomotopy([g], [f])

solve(H, [[-2.1 + 4.0im], [4.3im]])
```

and if you don't have a start system at hand we can use a total degree homotopy
```julia
using HomotopyContinuation
using DynamicPolynomials

@polyvar x
f = (x - 2.0) * (x - (2.5+ 4.0im))
H, startsolutions = totaldegree(StraightLineHomotopy, [f])

solve(H, startsolutions)
```
