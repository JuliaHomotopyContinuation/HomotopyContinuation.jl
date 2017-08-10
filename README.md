# HomotopyContinuation.jl [WIP]
[![Build Status](https://travis-ci.org/saschatimme/HomotopyContinuation.jl.svg?branch=master)](https://travis-ci.org/saschatimme/HomotopyContinuation.jl)
[![codecov](https://codecov.io/gh/saschatimme/HomotopyContinuation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/saschatimme/HomotopyContinuation.jl)

This is a package for solving systems of polynomials via homotopy continuation.

Getting started
-----------
Since this package is pre-release you have to install it via
```sh
Pkg.clone("git:://github.com/saschatimme/HomotopyContinuation.jl.git")
```
And here is a simple example:
```julia
using HomotopyContinuation
using TypedPolynomials

@polyvar x
F = [(x - 2.0) * (x - (2.5+ 4.0im))]
G = [(x - 4.3im) * (x + (2.1 - 4im))]

H = StraightLineHomotopy(G, F)

res = solve(H, [[-2.1 + 4.0im], [4.3im]], PredictorCorrector.Spherical())
```
