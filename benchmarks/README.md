# Benchmarks

This directory contains benchmarks for numerical kernels and polynomial systems.
The larger problems are too expensive for the normal test suite and should be
checked manually after major changes to the package.

Run comparisons sequentially on an idle machine, using the same Julia version
and dependency versions for each checkout. Check accuracy and success counts
alongside timings. The scripts below exclude compilation from their measurements.

## Numerical kernels

The following scripts measure backend kernels without changing solver options:

```sh
julia --project=. benchmarks/linear_algebra_kernels.jl
julia --project=. benchmarks/lu_factorization.jl
julia --project=. benchmarks/double_double_kernels.jl
julia --project=. benchmarks/modelkit_extended.jl
julia --project=. benchmarks/modelkit_extended.jl --unfused
```

The linear algebra and arithmetic scripts compare optimized operations with scalar or ordinary-arithmetic
baselines. `modelkit_extended.jl` measures double-double evaluation, Jacobians, and
third-order Taylor evaluation on four existing test systems. Its `--unfused` mode
disables the compound arithmetic specializations only in the benchmark process,
so both runs use the same precision. Add `--compiled` to compare generated
evaluators instead of interpreted ones. Kernel timings do not measure complete
tracking runs or batch solves.

## Path tracking

```sh
julia --project=. benchmarks/path_tracker.jl
```

`path_tracker.jl` measures complete tracking runs on regular paths, a system
requiring an extended-precision Jacobian, and the difficult Steiner fixture.
It reports minimum and median times after warm-up, endpoint errors, residuals,
and accepted/rejected step counts. It fails if a path does not reach its expected
endpoint.

To compare checkouts, run the same script with `--project=/path/to/checkout` for
each version. These measurements cover individual paths, not full batch solves
or certification.
