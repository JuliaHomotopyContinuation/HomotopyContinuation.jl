# Path tracking harness

Tracks fixed parameter-homotopy scenarios with the checkout's `Tracker` and `solve`, and
records robustness (failures, path jumping, certified solutions) and cost (steps, time)
per system, so that every change to tracking can be measured against a fixed reference.
The systems come from the solver-independent collection in `../systems` (see its README).

## Setup

```
julia --project=benchmarks/harness -e 'using Pkg; Pkg.instantiate()'
```

The environment uses the HomotopyContinuation.jl of this checkout (`[sources]` in
`Project.toml`, Julia ≥ 1.11).

Build the scenario caches once per machine (`benchmarks/cache/`, not committed):

```
julia --project=benchmarks/harness benchmarks/harness/reference.jl
```

The start data of every scenario (p₀, p₁ and all solutions at p₀) is the artifact
`benchmark_starts` (`Artifacts.toml`, format in `starts.jl`), downloaded on first use, so
every machine tracks the same paths in the same order. `reference.jl` builds
`benchmarks/cache/<scenario>.jls` from it: the start data and a reference run of every path
(compile = :all, single thread) with the HC.jl that built the cache (`hc_version`).
Existing caches are not rebuilt, so the reference stays fixed across commits; its timings
are only comparable on the machine that built it.

### Changing the start data

To add a system or scenario, add it to `load_system.jl` / `reference.jl` and run
`reference.jl`: scenarios that are not in the artifact get new start data (seeded random
parameters, start solutions from `solve`, or for hard systems from monodromy with
`solution_sets.jl`, run first). Then write a new artifact from the caches and bind it:

```
julia --project=benchmarks/harness benchmarks/harness/export_starts.jl <url> <tarball>
```

Upload the tarball so it downloads from `<url>` (a new asset name per version, e.g. a
GitHub release asset `benchmark_starts-v2.tar.gz`), check it with `certify.jl --starts`,
and commit `Artifacts.toml`. An artifact is never edited in place.

## Workflow

1. Commit the change.
2. Track all scenarios on the clean commit:
   `julia --project=benchmarks/harness benchmarks/harness/run.jl <label> [interpreted|compiled]`
   and, for scaling invariance, again with a third argument `scaled`: every system is
   tracked as G(y) = c .* F(D y) with fixed random cᵢ, dⱼ ∈ 10^U(-3, 3), applied around
   the unchanged program (inputs times d, outputs times c); rows get the system suffix
   `~scaled`, files the suffix `_scaled`. Endpoints are mapped back to the original
   coordinates, so certification and the comparison with the reference work unchanged.
   Trailing `option=value` arguments set `TrackerOptions` fields.
3. Certify the endpoints (the run's and the reference's) with interval arithmetic:
   `julia --threads=auto --project=benchmarks/harness benchmarks/harness/certify.jl <commit>_<label>_<mode>`
4. Add an entry to `../log/CHANGELOG.md` and commit `../log/runs.tsv`,
   `../log/certification.tsv` and the entry.

Labels must not contain `_`. Runs from a worktree with uncommitted changes under `src/` or
`benchmarks/harness/` get the commit `<sha>-dirty`; don't commit those rows.
`HC_BENCH_ONLY=<regex>` restricts a run to matching scenarios (smoke tests), and
`HC_BENCH_RUNS=<path>` / `HC_BENCH_CERTIFICATION=<path>` write the rows elsewhere.

## Scripts

- `run.jl`: tracks every scenario with `Tracker` (interpreted or compiled, optionally
  scaled) and compares each path with the reference.
- `solve.jl <label> [interpreted|compiled] [scaled]`: runs `solve` (tracking with the
  endgame) from the start solutions; rows in `../log/solve_runs.tsv`.
- `certify.jl`: certifies a run's endpoints (`certify.jl --starts`: the artifact's start
  sets, at their parameters).
- `certify_time.jl <run> [system ...]`: single-threaded time of `certify` per endpoint.
- `compare.jl <sha> <label-prefix> <reference-run>`: compares runs path by path with
  another run (e.g. of the commit a change starts from): paths that end at the same
  endpoint, fail, or jump visibly (duplicate endpoint), silently (another reference
  endpoint, no duplicate) or to an unknown point.
- `monodromy.jl <scenario> [loops] [seed]`: monodromy loops p₀ → p_a → p_b → p₀ through
  `solve` from a full solution set; each leg must return all solutions, certified and
  distinct, and each loop must permute the set.
- `load_system.jl`: stdlib-only loader for `../systems/*.toml` and the list of harness
  systems; `scaled.jl`: the `scaled` mode.
- `starts.jl`: reads the artifact `benchmark_starts`; `export_starts.jl` writes it;
  `solution_sets.jl`: monodromy solution sets for new start data of hard systems.

## Output

- `../log/runs.tsv`: one row per system per run, appended. Columns:
  - identity: `date`, `commit`, `label`, `mode` (interpreted/compiled), `julia`, `cpu`,
    `threads`, `system`, `paths`
  - robustness: `success`, `step_size_too_small`, `max_steps`, `accuracy_limit`,
    `ill_conditioned`, `other_failures` (failure counts by return code), `duplicates`
    (successful endpoints equal to an earlier one: a sign of path jumping),
    `mismatch_ref` (paths both runs finished but with different endpoints), `ref_success`
  - cost: `accepted`, `rejected` (steps, totals over all paths), `extended_precision`
    (paths that used extended precision), `total_s`, `median_us`, `p90_us` (wall time
    per path), `ref_accepted`, `ref_median_us`
  - `results_file`: per-path file in `../results/`
- `../results/<commit>_<label>_<mode>.tsv` (not committed): one row per path (return
  code, steps, time, final `t`, `accuracy`), with the reference's return code, steps and
  time for the same path alongside.
- `../log/certification.tsv`: one row per run, system and tracker (`run`, `reference`, or
  `start set`): `success` endpoints, how many are `certified`
  solutions, `distinct_certified` (the robust count of solutions found),
  `real_certified`, `distinct_real_certified`.
- `../cache/endpoints/<commit>_<label>_<mode>.jls` (not committed): the endpoints of a
  run, read by `certify.jl` and `compare.jl`.

## Scenarios

`HARNESS_SYSTEMS` in `load_system.jl` lists the systems; each is tracked between two
seeded random complex parameters. The hard systems (four-bar 8652, Fano quintic 2875,
Steiner 3264, tritangents 720) start from full solution sets computed with monodromy
(`solution_sets.jl`); each set is certified complete and distinct (`certify.jl --starts`).
Each hard system also has `<name>_g1`..`_g3`: the same start set tracked to three more
random complex targets (p₁ seeded 1, 2, 3). `steiner_conics_real` tracks the generic
Steiner set to the real five-conic configuration of HC.jl's certification test, which has
3264 real solutions; `steiner_conics_real_g<k>` tracks from the target of
`steiner_conics_g<k>` (the reference's endpoints there, all successful and distinct) to
that configuration. The systems from HC.jl issues also get `<name>_instance`: generic →
the stored instance of the issue.
