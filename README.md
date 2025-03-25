<img src="https://www.juliahomotopycontinuation.org/images/logo_transparent_bg.png" width="320px">

[![][docs-stable-img]][docs-stable-url] [![Run tests](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/actions/workflows/run_tests.yml/badge.svg)](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/actions/workflows/run_tests.yml)


**HomotopyContinuation.jl** is a Julia package for solving systems of polynomial equations by numerical homotopy continuation.

---

### **See [juliahomotopycontinuation.org](https://www.juliahomotopycontinuation.org) for installation instructions, full content overview and detailed documentation.**

---

## Basic usage

HomotopyContinuation.jl aims at having easy-to-understand top-level commands. Here is a simple example:

```julia
using HomotopyContinuation
@var x y; # declare the variables x and y
F = System([x^2+2y, y^2-2])
result = solve(F)
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
  title={{H}omotopy{C}ontinuation.jl: {A} {P}ackage for {H}omotopy {C}ontinuation in {J}ulia},
  author={Breiding, Paul and Timme, Sascha},
  booktitle={International Congress on Mathematical Software},
  pages={458--465},
  year={2018},
  organization={Springer}
}
```

A preprint of this paper is [freely available](https://arxiv.org/abs/1711.10911).

[docs-stable-img]: https://img.shields.io/badge/docs-online-blue.svg
[docs-stable-url]: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable
