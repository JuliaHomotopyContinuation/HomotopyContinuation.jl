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

- `linear_algebra_kernels.jl` compares the paired-column triangular solves
  `ldiv_upper!` and `ldiv_unit_lower!` with a scalar-column baseline. They are
  used by the condition estimates and the QR solves of overdetermined systems.
- `lu_factorization.jl` times the LU factorization and solve of square
  `MatrixWorkspace`s, which hold the Jacobians during path tracking and store the
  real and imaginary parts of their factors separately.
- `double_double_kernels.jl` compares ModelKit's fused double-double operations
with ordinary complex double-double arithmetic.
- `modelkit_extended.jl` measures double-double evaluation, Jacobians, and
third-order Taylor evaluation on four existing test systems. Its `--unfused`
mode disables the compound arithmetic specializations only in the benchmark
process, so both runs use the same precision. Add `--compiled` to measure
`CompiledSystem`s instead of interpreted ones. Compiled systems evaluate
double-double Jacobians and Taylor series with an interpreter, so in this mode
only the evaluation column measures generated code.


## Path tracking

```sh
julia --project=. benchmarks/path_tracker.jl
```

`path_tracker.jl` measures complete tracking runs on regular paths, a system
requiring an extended-precision Jacobian, and the difficult Steiner fixture.
It reports minimum and median times after warm-up, endpoint errors, residuals,
and accepted/rejected step counts. It fails if a path does not reach its expected
endpoint.


## Other scripts

- `run_system_benchmark.jl` and `judge_system_benchmark.jl` define
  `run_benchmark(file)` and `judge_results(file_pr, file_main)`. They use
  BenchmarkTools to time `evaluate!`, `evaluate_and_jacobian!`, and `taylor!` of
  the systems in `test/test_systems.jl`, interpreted and compiled, with the
  tuned parameters in `params.json`, and compare two saved runs. BenchmarkTools
  is not a dependency of this package and must be installed separately.
- `cyclooctane.jl`, `tritangents.jl`, and `bio-chemical-rection-networks.jl` are
  older problem scripts written with DynamicPolynomials. `tritangents.jl` also
  needs PolynomialTestSystems, which is not a dependency of this package.
- `generate_path_tracker_seed_data.jl` includes `path_tracker_benchmark_utils.jl`,
  which is not in the repository, so it does not run.
