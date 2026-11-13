# Benchmark polynomial systems

A collection of polynomial systems for benchmarking **parameter-homotopy path
tracking** that does not depend on any one solver. Each system comes with its provenance and its expected number of
solutions. Every system is one TOML file that HomotopyContinuation.jl (HC.jl) or any other solver
can read. The equations use plain infix syntax that Julia's `Meta.parse` accepts.

```
benchmarks/systems/
  tier_a/   realistic parameter homotopies (kinematics, vision, chemistry, power flow, geometry)
  tier_b/   stress / scaling families (coefficient-parameter families, path-jumping constructions)
  tier_c/   endgame targets: non-generic target parameters with singular or infinite endpoints
  tier_d/   large systems: many unknowns, many terms, many parameters (generated, not committed:
            run tier_d/generate.jl)
  hcjl/     every system from HC.jl's test suite and benchmarks/ that is not already in tier_a-c
  validate.jl
```

Check the collection with (Julia standard library only):

```
julia benchmarks/systems/validate.jl
```

It checks every file against the schema below. Equations must parse, use only the allowed operations, and use
only declared names. Systems must be square unless marked otherwise. Value vectors must have the right lengths.
Every equation must evaluate to a finite number at a random complex point.

## File format

```toml
name = "stewart_gough_fk"          # == file name
tier = "A"                         # "A" | "B" | "C" | "D" | "HC" (== directory)
description = "..."
tags = ["kinematics", ...]
variables = ["n1", "n2", ...]      # unknowns, in order
parameters = ["A21", ...]          # may be empty (see [parametrization])
equations = ["n1^2 + n2^2 + n3^2 - 1", ...]

# optional top-level keys
variable_groups = [[...], [...]]   # multi-homogeneous structure (variables are listed in group order)
square = false                     # only for over/underdetermined systems
example_parameters = ["0.645275", ...]      # a realistic instance from the source
start_parameters   = ["0.1+2.0im", ...]     # start pair from the source (tier C / hcjl test cases)
start_solution     = [...]
target_parameters  = [...]                  # target instance (required in tier C)
target_solution    = [...]
generator = "..."                  # how to produce other members of a family (e.g. katsura-n)
benchmark_paths = "..."            # which parameter paths to track (path-jumping constructions)

[source]
citation = "..."                   # paper(s) for the system and its count
url = "..."                        # pinned link to the file the equations were taken from
file = "Repo@commit:path:lines"    # same, as text

[solutions]
generic_nonsingular = 40           # nonsingular solutions for generic (complex) parameters
generic_nonsingular_formula = "2^N = 2^100"  # instead, when the count does not fit in a TOML integer
count_source = "..."               # where that number comes from
target_endpoints = "..."           # tier C: what happens at target_parameters
expected_outcome = "..."           # single-path / path-jumping cases: what a correct tracker returns
zero_dimensional = false           # only for systems where no count applies
count_known = false                # count not known (hcjl/rigid_multiview)

[parametrization]
kind = "explicit"                  # "explicit" | "coefficients" | "pointer" | "fixed"
description = "..."
generic_parameters = [...]         # kind = "coefficients": optional explicit generic family
generic_equations = [...]

[verification]                     # what was checked in this collection (see below)
status = "verified"                # "verified" | "partial" | "not-verified" | "not-needed"
method = "..."
hc_version = "2.22.4"
verified_count = 40

[hcjl]                             # hcjl/ files: how HC.jl itself parametrizes and solves the system
[start_data]                       # pointers to GPU-HC start/target data files
[external]                         # kind = "pointer": where the data lives
```

**Equation syntax.** Use `+ - * / ^` with nonnegative integer literal exponents. Numbers may be
integers, decimals (including `1e-3`), quotients `29/16`, and complex numbers written with `im`
(`0.5-1.2im`). `sqrt` appears only in parameter expressions (`hcjl/sqrt_parameters`). Division
appears only in the rational HC.jl systems `hcjl/rigid_multiview` and `hcjl/small_rational`. Names
must be ASCII Julia identifiers and must not be `im`. To parse a file, read it with `TOML.parsefile`,
`Meta.parse` each equation, and map symbols to your own variable type (`validate.jl` has a 10-line
evaluator to copy).

**Parametrization kinds.**

- `explicit`: `parameters` is the parameter vector of the family; `equations` depend on it.
- `coefficients`: `parameters` is empty, and the stored equations are one member of the
  coefficient-parameter family. The family is defined by multiplying every monomial term of the
  expanded equations by an independent parameter. `generic_equations` / `generic_parameters` give
  that family explicitly when it is small. The generic count of the family is
  `generic_nonsingular`.
- `pointer`: the data is too large or is already well curated elsewhere. `[external]` names the
  repository, the revision and the files.
- `fixed`: not a parameter family (for example HC.jl's overdetermined `minors`).

**Parameter order matters.** Parameters are listed in the order of the source (for HC.jl systems,
`System(...).parameters` order). Value vectors such as `start_parameters` use that order.

**Verification.** HC.jl 2.22.4 (Julia 1.13.3) was used in a scratch environment. Equations taken
from HC.jl were compared with the original HC.jl functions at random complex points. Equations taken
from PHCpack were compared with the PHCpack Demo files. In both cases the agreement is about 1e-16
(stewgou40 agrees to 2e-10 because PHCpack rounds its coefficients). The counts were then
recomputed in one of three ways:

- `monodromy_solve`: a lower bound that stops after 10-20 loops without progress.
- a polyhedral or total-degree `solve`.
- for GPU-HC problems, monodromy from one GPU-HC start solution, plus a check that all GPU-HC start
  solutions have a small residual.

`[verification]` records the method used for each system.

## Systems

`#x` / `#p`: number of unknowns / parameters. `count`: generic nonsingular solution count.
`HC.jl ref` points to the HC.jl source the system was taken from, at commit `0cbf1e27` (HC.jl 2.22.4).

| tier | system | #x | #p | count | verified (HC.jl 2.22.4) | tags | source | HC.jl ref |
|---|---|---|---|---|---|---|---|---|
| A | [`alt_nine_point_four_bar`](tier_a/alt_nine_point_four_bar.toml) | 24 | 16 | 8652 | yes (8652) | kinematics, mechanism-synthesis, parameter-homotopy | HomotopyContinuation.jl:test/test_systems.jl:954-979 (four_bar) | test/test_systems.jl:954-979 |
| A | [`bacillus_subtilis`](tier_a/bacillus_subtilis.toml) | 10 | 1 | 44 | yes (44) | chemistry, reaction-network, steady-states | HomotopyContinuation.jl:test/test_systems.jl:79-129 (bacillus) | test/test_systems.jl:79-129 |
| A | [`chicago_trifocal`](tier_a/chicago_trifocal.toml) | 30 | 33 | 312 | yes (312) | vision, minimal-problem, trifocal | Homotopy-Continuation-Tracker-on-GPU:auto-data-gen-tools/p2c-problems/sys_trifocal_2op1p_30x30.m |  |
| A | [`five_point_relative_pose`](tier_a/five_point_relative_pose.toml) | 6 | 20 | 40 | yes (40) | vision, minimal-problem, relative-pose | Homotopy-Continuation-Tracker-on-GPU:auto-data-gen-tools/p2c-problems/sys_5pt_rel_pos_alg_form_quat.m |  |
| A | [`lossless_powerflow_k4`](tier_a/lossless_powerflow_k4.toml) | 6 | 9 | 20 | yes (20) | power-flow, kuramoto, parameter-homotopy | https://doi.org/10.1109/TCS.1982.1085093 |  |
| A | [`six_revolute_ik`](tier_a/six_revolute_ik.toml) | 12 | 15 | 16 | yes (16) | kinematics, robotics, inverse-kinematics | HomotopyContinuation.jl:test/test_systems.jl:695-728 (six_revolute) | test/test_systems.jl:695-728 |
| A | [`steiner_conics`](tier_a/steiner_conics.toml) | 15 | 30 | 3264 | yes (3264) | enumerative-geometry, parameter-homotopy, large-count | HomotopyContinuation.jl:test/test_systems.jl:873-900 (steiner) | test/test_systems.jl:873-900 |
| A | [`stewart_gough_fk`](tier_a/stewart_gough_fk.toml) | 9 | 29 | 40 | yes (40) | kinematics, robotics, forward-kinematics | PHCpack:src/Demo/stewgou40 |  |
| A | [`triangulation_3view`](tier_a/triangulation_3view.toml) | 8 | 24 | 31 | yes (31) | vision, triangulation, critical-points | Homotopy-Continuation-Tracker-on-GPU:auto-data-gen-tools/p2c-problems/sys_3view_triangulation.m |  |
| A | [`triangulation_4view`](tier_a/triangulation_4view.toml) | 11 | 35 | 142 | yes (142) | vision, triangulation, critical-points | Homotopy-Continuation-Tracker-on-GPU:auto-data-gen-tools/p2c-problems/sys_4view_triangulation.m |  |
| B | [`algpath_clustered_roots`](tier_b/algpath_clustered_roots.toml) | 1 | 1 | nc*cs (instances: 25, 15, 50, 100, 100) | n/a | path-jumping, clustered-roots, stress | algpath-bench:gen.sage (class ClusteredHomotopy) |  |
| B | [`coinciding_paths_no_branch_point`](tier_b/coinciding_paths_no_branch_point.toml) | 1 | 3 | 2 | n/a | path-jumping, stress, univariate | https://arxiv.org/abs/1909.04984 |  |
| B | [`cyclic7`](tier_b/cyclic7.toml) | 7 | 0 | 924 | yes (924) | stress, coefficient-parameter, symmetric | HomotopyContinuation.jl:test/test_systems.jl:72-77 (cyclic) | test/test_systems.jl:72-77 |
| B | [`dense_univariate_1000`](tier_b/dense_univariate_1000.toml) | 1 | 1001 | 1000 | n/a | stress, scaling, univariate | algpath-bench:data/linear/dense-1000/all/ |  |
| B | [`dense_univariate_500`](tier_b/dense_univariate_500.toml) | 1 | 501 | 500 | n/a | stress, scaling, univariate | algpath-bench:data/linear/dense-500/all/ |  |
| B | [`eco8`](tier_b/eco8.toml) | 8 | 0 | 64 | yes (64) | stress, coefficient-parameter, diverging-paths | PHCpack:src/Demo/eco8 |  |
| B | [`hyperbola_path_jumping`](tier_b/hyperbola_path_jumping.toml) | 1 | 2 | 2 | n/a | path-jumping, stress, univariate | https://arxiv.org/abs/1909.04984 |  |
| B | [`katsura10`](tier_b/katsura10.toml) | 11 | 0 | 1024 | yes (1024) | stress, scaling, coefficient-parameter | PolynomialTestSystems.jl:src/systems.jl (katsura) |  |
| B | [`katsura12`](tier_b/katsura12.toml) | 13 | 0 | 4096 | yes (4096) | stress, scaling, coefficient-parameter | PolynomialTestSystems.jl:src/systems.jl (katsura) |  |
| B | [`noon5`](tier_b/noon5.toml) | 5 | 0 | 233 | yes (233) | stress, coefficient-parameter, symmetric | PHCpack:src/Demo/noon5 |  |
| C | [`double_root_issue614`](tier_c/double_root_issue614.toml) | 1 | 2 | 2 | partial (2) | endgame, singular-endpoint, multiplicity | https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/614 |  |
| C | [`griewank_osborne`](tier_c/griewank_osborne.toml) | 2 | 9 | 3 | yes (3) | endgame, singular-endpoint, multiplicity | PolynomialTestSystems.jl:src/systems.jl (griewank_osborne) |  |
| C | [`hyperbolic_6_6`](tier_c/hyperbolic_6_6.toml) | 2 | 9 | 12 | partial (12) | endgame, singular-endpoint, multiplicity | HomotopyContinuation.jl:test/endgame_test.jl:7-22 | test/endgame_test.jl:7-22 |
| C | [`infinite_endpoints_quadrics`](tier_c/infinite_endpoints_quadrics.toml) | 2 | 12 | 4 | yes (4) | endgame, infinite-endpoint, diverging-paths | HomotopyContinuation.jl:test/endgame_test.jl:54-60 | test/endgame_test.jl:54-60 |
| HC | [`biochem_network1`](hcjl/biochem_network1.toml) | 3 | 10 | 3 | yes (3) | chemistry, reaction-network, steady-states | HomotopyContinuation.jl:benchmarks/bio-chemical-rection-networks.jl:19-44 | benchmarks/bio-chemical-rection-networks.jl:19-44 |
| HC | [`biochem_network2`](hcjl/biochem_network2.toml) | 3 | 8 | 6 | yes (6) | chemistry, reaction-network, steady-states | HomotopyContinuation.jl:benchmarks/bio-chemical-rection-networks.jl:46-104 | benchmarks/bio-chemical-rection-networks.jl:46-104 |
| HC | [`biochem_network4`](hcjl/biochem_network4.toml) | 1 | 3 | 3 | yes (3) | chemistry, reaction-network, univariate | HomotopyContinuation.jl:benchmarks/bio-chemical-rection-networks.jl:106-125 | benchmarks/bio-chemical-rection-networks.jl:106-125 |
| HC | [`cyclic5`](hcjl/cyclic5.toml) | 5 | 0 | 70 | yes (70) | coefficient-parameter, symmetric, polyhedral | HomotopyContinuation.jl:test/test_systems.jl:72-77 (cyclic(n), n = 5) | test/test_systems.jl:72-77 |
| HC | [`cyclooctane`](hcjl/cyclooctane.toml) | 17 | 36 | 1408 | partial | chemistry, molecular-conformation, witness-set | HomotopyContinuation.jl:test/test_systems.jl:57-70 (cyclo, 15 equations, no slice) | test/test_systems.jl:57-70; benchmarks/cyclooctane.jl:1-29 |
| HC | [`fano_quintic_lines`](hcjl/fano_quintic_lines.toml) | 6 | 125 | 2875 | yes (2875) | enumerative-geometry, parameter-homotopy, large-count | HomotopyContinuation.jl:test/test_systems.jl:730-870 (fano_quintic) | test/test_systems.jl:730-870 |
| HC | [`four_bar_tracker_case`](hcjl/four_bar_tracker_case.toml) | 24 | 16 | 8652 | yes (single path) | kinematics, single-path, path-jumping | HomotopyContinuation.jl:test/test_cases/four_bar.jl:1-97 | test/test_cases/four_bar.jl:1-97 |
| HC | [`minors_3x5`](hcjl/minors_3x5.toml) | 3 | 0 | 80 | yes (80) | overdetermined, fixed, excess-solutions | HomotopyContinuation.jl:test/test_systems.jl:131-692 (minors) | test/test_systems.jl:131-692 |
| HC | [`moments3_gaussian_mixture`](hcjl/moments3_gaussian_mixture.toml) | 9 | 9 | 1350 | yes (1350) | statistics, parameter-homotopy, symmetric | HomotopyContinuation.jl:test/test_systems.jl:903-952 (moments3) | test/test_systems.jl:903-952 |
| HC | [`rigid_multiview`](hcjl/rigid_multiview.toml) | 7 | 29 | n/a | no | vision, triangulation, critical-points | HomotopyContinuation.jl:test/test_systems.jl:1-29 (RigidMultiView | test/test_systems.jl:1-29 |
| HC | [`small_rational`](hcjl/small_rational.toml) | 6 | 8 | n/a | n/a | rational, evaluation-only, not-square | HomotopyContinuation.jl:test/test_systems.jl:981-995 (small_rational) | test/test_systems.jl:981-995 |
| HC | [`sqrt_parameters`](hcjl/sqrt_parameters.toml) | 2 | 2 | 6 | yes (6) | non-polynomial-parameters, sqrt, small | HomotopyContinuation.jl:test/test_systems.jl:997-1001 (sqrt_parameters) | test/test_systems.jl:997-1001 |
| HC | [`steiner_higher_prec_case`](hcjl/steiner_higher_prec_case.toml) | 15 | 30 | 3264 | yes (single path) | enumerative-geometry, single-path, near-singular | HomotopyContinuation.jl:test/test_cases/steiner_higher_prec.jl:1-153 | test/test_cases/steiner_higher_prec.jl:1-153 |
| HC | [`tritangents`](hcjl/tritangents.toml) | 12 | 20 | 720 | yes (720) | enumerative-geometry, parameter-homotopy, symmetric | HomotopyContinuation.jl:test/test_systems.jl:31-53 (Tritangtens) | test/test_systems.jl:31-53 |

### Tier D: large systems

Tier D stresses evaluation and tracking at scale: many unknowns (up to 100), many terms (up to
11440), many parameters (up to 4004), or very large solution counts. The files are not
committed: [`tier_d/generate.jl`](tier_d/generate.jl) writes all of them, byte for byte the
same on every run (standard library only; run `julia benchmarks/systems/tier_d/generate.jl`). Its functions take the size as an argument, so other
members of each family (katsura-n, cyclic-n, eco-n, noon-n, dense n/d, K_n power flow, BVP with N
nodes) can be produced the same way. The generator reproduces `tier_b/katsura12`, `tier_b/cyclic7`,
`tier_b/eco8`, `tier_b/noon5` and `tier_a/lossless_powerflow_k4` exactly. Every file is fully
expanded (the largest, `fano_septic_lines`, is 250 KB), so `validate.jl` checks all of them like the
other tiers.

`#terms`: number of monomials in the unknowns, summed over all equations, after expansion
(parameters count as coefficients).

| system | #x | #p | #terms | max deg | count | verified (HC.jl 2.22.4) | tags | source |
|---|---|---|---|---|---|---|---|---|
| `katsura15` | 16 | 0 | 216 | 2 | 32768 | yes (32768) | sparse, coefficient-parameter | PolynomialTestSystems.jl:src/systems.jl (katsura) |
| `katsura20` | 21 | 0 | 362 | 2 | 1048576 | no | sparse, coefficient-parameter | PolynomialTestSystems.jl:src/systems.jl (katsura) |
| `cyclic10` | 10 | 0 | 92 | 10 | 35940 (family) | yes (35940) | symmetric, polyhedral, high-degree | PHCpack:src/Demo/cyclic10 |
| `eco15` | 15 | 0 | 134 | 3 | 8192 | yes (8192) | diverging-paths, coefficient-parameter | PHCpack:src/Demo/eco8 (family) |
| `noon9` | 9 | 0 | 90 | 3 | 19665 | yes (19665) | symmetric, coefficient-parameter | PHCpack:src/Demo/noon5 (family) |
| `bvp_quadratic_n100` | 100 | 101 | 498 | 2 | 2^100 | partial (N = 8, 10) | boundary-value, tridiagonal, single-path | https://doi.org/10.1007/s00607-005-0132-4 |
| `dense_n6_d4` | 6 | 1260 | 1260 | 4 | 4096 | n/a (Bezout) | dense, many-terms, many-parameters | Bezout's theorem |
| `dense_n8_d4` | 8 | 3960 | 3960 | 4 | 65536 | n/a (Bezout) | dense, many-terms, many-parameters | Bezout's theorem |
| `dense_n4_d10` | 4 | 4004 | 4004 | 10 | 10000 | n/a (Bezout) | dense, many-terms, high-degree | Bezout's theorem |
| `lossless_powerflow_k10` | 18 | 54 | 189 | 2 | 48620 | yes (48620) | power-flow, kuramoto, parameter-homotopy | https://doi.org/10.1109/TCS.1982.1085093 |
| `fano_septic_lines` | 8 | 791 | 11440 | 7 | 698005 | no | enumerative-geometry, many-terms, many-parameters | https://oeis.org/A027363 |

Notes on tier D:

- **Counts.** For `cyclic10` the count is for the coefficient family (mixed volume 35940); the
  system with all coefficients 1 has 34940 isolated solutions (Dai-Kim-Kojima 2003).
  `katsura15` and `katsura20` have torus mixed volume below 2^n (32730 for n = 15), like
  `katsura12`. `eco15` and `noon9`: the counts 2^(n-2) and 3^n - 2n match the HC.jl mixed volume for
  every n = 8..16 and n = 5..9.
- **`bvp_quadratic_n100`** is meant for tracking a few paths (for example the real solutions of the
  unforced problem at `example_parameters`) in a large sparse system. Its count 2^100 does not fit
  in a TOML integer and is stored as `generic_nonsingular_formula`. The forcing parameters `f1..f100`
  make the count provable; see the file.
- **`lossless_powerflow_k10`**: the mixed volume is 257536, more than five times the 48620
  solutions, so polyhedral homotopies waste most paths (as for `tier_a/lossless_powerflow_k4`:
  40 vs. 20). Use monodromy or a parameter homotopy.
- **Dense systems** are the coefficient-parameter families themselves; the total-degree start
  system is a member (see each file).
- **Considered but not included**: the Chandrasekhar H-equation for n = 16 (PHCpack `chandra4-6`
  family, mixed volume 2^(n-1)). HC.jl 2.22.4 polyhedral finds only 124 of 128 solutions at n = 8,
  823 of 2048 at n = 12 and about 3000 of 32768 at n = 16 (solutions with norm 1e5 and larger), and
  monodromy did not get past one solution, so the count could not be checked.

### Notes on specific systems

- **Monodromy undercounts reducible problems**: the solutions of `hcjl/cyclooctane` lie on a
  reducible variety (monodromy from one point finds 632 of 1408). `hcjl/biochem_network4` has
  `x = 0` as a separate component (monodromy finds 2 of 3). Count these systems with a
  complete start system.
- **Endgame findings with HC.jl 2.22.4** (tier C): the parameter homotopy from a generic member
  of `hyperbolic_6_6` to its target fails on all 12 paths, although HC.jl's total-degree homotopy
  handles the target. Issue #614 (`double_root_issue614`) still occasionally fails, in 2 of 10
  runs.
- **bacillus**: one polyhedral run at a random complex `T` returned 43 instead of 44, so this
  badly scaled system is a useful robustness test.
- **Alt's four-bar**: `tier_a/alt_nine_point_four_bar` and `hcjl/four_bar_tracker_case` contain
  the same 24x24 system (HC.jl `four_bar()`). The tier A file is the full 8652-path parameter
  homotopy; 1442 is the count up to the 6-fold cognate/relabeling symmetry. The hcjl file is HC.jl's
  single hard path from complex parameters `p` to physical parameters `q`
  (`test/test_cases/four_bar.jl`), with its start solution. It is stored with 8 digits only, so
  its residual is 4.4e-7: refine it with Newton's method before tracking.
- **Steiner**: `hcjl/steiner_higher_prec_case` is HC.jl's single path (complex `p` -> real `q`)
  that passes close to a singularity and needs extended-precision residuals for 1-2 steps.
- **GPU-HC problems** (5-point, Chicago, triangulation): the equations are copied from the
  MATLAB generator files. GPU-HC also ships start parameters and start solutions (`[start_data]`,
  paths relative to `GPU-HC/problems/`). The triangulation counts (31, 142) are for GPU-HC's
  formulation with consecutive-pair epipolar constraints; exact 3-view triangulation has 47
  stationary points.
- **katsura12**: HC.jl's `mixed_volume` reports 4020 (torus count). The affine count is 4096,
  because 76 solutions have a zero coordinate. A purely toric polyhedral start system therefore
  misses solutions.
- **Pointers only**: `tier_b/algpath_clustered_roots` (algpath-bench JSON instances). The dense
  univariate degree 500/1000 families are stored explicitly as generic forms, and the concrete
  algpath instances are referenced.
- **Not included**: no reliable copyable source with known counts was found for larger power-flow
  instances (for example IEEE test cases). `tier_a/lossless_powerflow_k4` and
  `tier_d/lossless_powerflow_k10` use a complete-graph formulation written for this collection.
  The larger GPU-HC/MINUS vision problems (e.g. `generalized_3views_6lines`, 600 solutions; `uncalibrated_trifocal_rel_pos_CH`,
  5344 tracks) are not copied; see `auto-data-gen-tools/p2c-problems/` in the GPU-HC repository.
