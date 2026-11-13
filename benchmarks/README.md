# Benchmarks

This file contains different scripts of polynomial systems and problems we encountered
during the development of this package.
These problems take too much time to include in our normal test suite, but nonetheless
these should be checked manually after major changes in the package.

- `systems/`: a solver-independent collection of benchmark polynomial systems, one TOML
  file per system with provenance and solution count (see `systems/README.md`).
- `harness/`: the path tracking harness, which tracks fixed scenarios built from those
  systems and logs robustness and cost per commit to `log/` (see `harness/README.md`).
- `log/`: the harness results of logged commits (`CHANGELOG.md` summarizes them).

## Running the path tracking harness

The harness tracks 40 fixed scenarios (80,875 paths, from 2 to 8652 per scenario) with
the `Tracker` of this checkout and compares every path with a fixed reference run. All
commands run from the repository root; Julia ≥ 1.11.

**Once per machine** (about 10 minutes): install the environment and build the scenario
caches in `benchmarks/cache/` (not committed). The start solutions are downloaded as the
artifact `benchmark_starts`; the reference run of every path is tracked locally, so its
timings belong to this machine.

```sh
julia --project=benchmarks/harness -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/harness benchmarks/harness/reference.jl
```

**Per change**: commit, then track all scenarios and certify the endpoints. `<label>` is
any name without `_`.

```sh
julia --project=benchmarks/harness benchmarks/harness/run.jl <label> compiled
julia --project=benchmarks/harness benchmarks/harness/run.jl <label> interpreted
julia --project=benchmarks/harness benchmarks/harness/run.jl <label> interpreted scaled
julia --threads=auto --project=benchmarks/harness benchmarks/harness/certify.jl <commit>_<label>_compiled
```

`run.jl` prints one line per scenario (paths, successes, duplicate endpoints, endpoints
that differ from the reference, steps and median time per path, next to the reference's)
and appends the same to `log/runs.tsv`; `certify.jl` appends the number of certified
distinct solutions to `log/certification.tsv`. Compare with the baseline rows in those
files (or run `compare.jl` against the baseline run, see `harness/README.md`), add an
entry to `log/CHANGELOG.md`, and commit the log.

For a quick check while iterating, restrict the scenarios with `HC_BENCH_ONLY` and write
the rows somewhere else (runs with uncommitted changes in `src/` get the commit
`<sha>-dirty` and do not belong in the log):

```sh
HC_BENCH_ONLY='^(katsura10|steiner_conics)\.' HC_BENCH_RUNS=/tmp/runs.tsv \
  julia --project=benchmarks/harness benchmarks/harness/run.jl wip compiled
```

`harness/README.md` describes all scripts (`solve.jl`, `compare.jl`, `monodromy.jl`, ...),
the output columns, the scenarios, and how to change the start data.
