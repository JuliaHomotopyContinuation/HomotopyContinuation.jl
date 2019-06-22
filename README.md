<img src="https://www.juliahomotopycontinuation.org/images/logo_transparent_bg.png" width="320px">

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-dev-img]][docs-dev-url] | [![Codecov branch][codecov-img]][codecov-url]|

**HomotopyContinuation.jl** is a Julia package for solving systems of polynomial equations by numerical homotopy continuation.

---

### **See [juliahomotopycontinuation.org](https://www.juliahomotopycontinuation.org) for installation instructions, full content overview and detailed documentation.**

---

## Basic usage

HomotopyContinuation.jl aims at having easy-to-understand top-level commands. Here is a simple example:

```julia
using HomotopyContinuation
@polyvar x y; # declare the variables x and y
result = solve([x^2+2y, y^2-2])
```
```
Result with 4 solutions
==================================
• 4 non-singular solutions (2 real)
• 0 singular solutions (0 real)
• 4 paths tracked
• random seed: 902575
```

For more see [our user guides](https://www.juliahomotopycontinuation.org/guides/).


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
