<img src="https://i.imgur.com/8ycOn14.png" width="320px">

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-dev-img]][docs-dev-url] | [![Codecov branch][codecov-img]][codecov-url]|

HomotopyContinuation.jl is available through the Julia package manager by


```julia-repl
pkg> add HomotopyContinuation
```


you can enter the Julia package manager by pressing `]` in the REPL.


## Basic usage

HomotopyContinuation.jl aims at having easy-to-understand top-level commands. For instance, suppose we want to solve the following polynomial system

```math
f= [x^2+y, y^2-1].  
```


This can be accomplished as follows


```julia
using HomotopyContinuation
@polyvar x y; # declare the variables x and y
result = solve([x^2+2y, y^2-2])
```

```
AffineResult with 4 tracked paths
==================================
• 4 non-singular finite solutions (2 real)
• 0 singular finite solutions (0 real)
• 0 solutions at infinity
• 0 failed paths
• random seed: 751138
```


Let us see what is the information that we get. Four paths were attempted to be solved, four of which were completed successfully. Since we tried to solve an affine system, the algorithm checks whether there are solutions at infinity: in this case there are none. With *solutions at infinity* we mean solutions of the [homogenization](https://en.wikipedia.org/wiki/Homogeneous_polynomial#Homogenization) of the system which are no solutions of the affine system. None of the solutions are singular and two of them are real. To access the first solution in the array we write


```julia
# Assume we are only interested in the *real* solutions. We can obtain these by
real(result)
```

```
2-element Array{HomotopyContinuation.PathResult{Complex{Float64},Float64,Complex{Float64}},1}:
 • returncode: success
 • solution: Complex{Float64}[1.68179+5.55112e-17im, -1.41421-1.38778e-17im]
 • residual: 1.130e-16
 • pathnumber: 3

 • returncode: success
 • solution: Complex{Float64}[-1.68179-2.77556e-17im, -1.41421+0.0im]
 • residual: 1.119e-16
 • pathnumber: 4
```


where we can see that there are 2 real solutions, `(⁴√(2³),-√2)` and `(-⁴√(2³), -√2)`. We also can look into more detail into the first result by


```julia
real(result)[1]
```

```
PathResult
==========
 • returncode: success
 • solution: Complex{Float64}[1.68179+5.55112e-17im, -1.41421-1.38778e-17im]
 • residual: 1.130e-16
 • condition_number: 1.640e+00
 • windingnumber: 1

 • pathnumber: 3
 • start_solution: Complex{Float64}[1.0+0.0im, -1.0+1.22465e-16im]

 • t: 0.0
 • iterations: 6
 • npredictions: 3
```


The returncode tells us that the pathtracking was successfull. What do the other entries of that table tell us? Let us consider the most relevant:


  * `solution`: the zero that is computed (here it is `(⁴√(2³),-√2)`).
  * `singular`: boolean that shows whether the zero is singular.
  * `residual`: the computed value of ``|f([1.68179+1.27738e-17im, -1.41421-1.18454e-17im])|``.
  * `windingnumber`: This is a lower bound on the multiplicity of the solution. A windingnumber greater than 1 indicates that the solution is singular.


To extract the solutions you can do


```julia
solutions(result)
```

```
4-element Array{Array{Complex{Float64},1},1}:
 [1.11022e-16+1.68179im, 1.41421+0.0im]
 [-1.11022e-16-1.68179im, 1.41421+1.11022e-16im]
 [1.68179+5.55112e-17im, -1.41421-1.38778e-17im]
 [-1.68179-2.77556e-17im, -1.41421+0.0im]
```


or if you only want real solutions


```julia
realsolutions(result)
```

```
2-element Array{Array{Float64,1},1}:
 [1.68179, -1.41421]
 [-1.68179, -1.41421]
```



## Citing HomotopyContinuation.jl
If you find HomotopyContinuation.jl useful in your work, we kindly request that you cite the [following paper](https://link.springer.com/chapter/10.1007/978-3-319-96418-8_54):

```latex
@inproceedings{HomotopyContinuation.jl,
  title={HomotopyContinuation.jl: A Package for Homotopy Continuation in Julia},
  author={Breiding, Paul and Timme, Sascha},
  booktitle={International Congress on Mathematical Software},
  pages={458--465},
  year={2018},
  organization={Springer}
}
```

A preprint of this paper is [freely available](https://arxiv.org/abs/1711.10911).

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable
[docs-dev-url]: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/dev

[build-img]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl
[codecov-img]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl
