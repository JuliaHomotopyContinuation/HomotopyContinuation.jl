# Benchmarks

This file contains different scripts of polynomial systems and problems we encountered
during the development of this package.
These problems take too much time to include in our normal test suite, but nonetheless
these should be checked manually after major changes in the package.

The following scripts measure backend kernels without changing solver options:

```sh
julia --project=. benchmarks/linear_algebra_kernels.jl
julia --project=. benchmarks/double_double_kernels.jl
julia --project=. benchmarks/modelkit_extended.jl
julia --project=. benchmarks/modelkit_extended.jl --unfused
```

The first two compare optimized operations with scalar or ordinary-arithmetic
baselines. `modelkit_extended.jl` measures double-double evaluation, Jacobians, and
third-order Taylor evaluation on four existing test systems. Its `--unfused` mode
disables the compound arithmetic specializations only in the benchmark process,
so both runs use the same precision. Add `--compiled` to compare generated
evaluators instead of interpreted ones. Run comparisons sequentially on an idle
machine; the reported warmed kernel times exclude compilation and are not
end-to-end solver speedups.
