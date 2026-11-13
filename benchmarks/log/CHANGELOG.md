# Harness log

What each logged commit changed in tracking and what it did to the harness numbers
(`runs.tsv`, `certification.tsv`, `solve_runs.tsv`; how to run: `../README.md`). Newest
first. Totals are over all 40 scenarios (80,875 paths); "certified" is the number of
distinct certified solutions over all scenarios (`certify.jl`).

## 40de5a31 baseline: HC.jl 2.23, before any harness-driven change

Apple M4 Max, Julia 1.13.3, one thread. Reference: the caches built with HC.jl 2.23
(`compile = :all`), so the compiled run repeats the reference exactly (same successes,
endpoints and step counts), which shows that runs are deterministic.

| run | success | step size too small | max steps | accuracy limit | duplicates | ≠ reference | accepted steps | tracking time | certified |
|---|---|---|---|---|---|---|---|---|---|
| compiled | 80,515 | 104 | 169 | 87 | 1 | 0 | 13,216,898 | 193.1 s | 80,296 |
| interpreted | 80,523 | 43 | 211 | 98 | 2 | 142 | 13,705,567 | 239.7 s | 80,297 |
| interpreted, scaled | 80,330 | 104 | 357 | 84 | 1 | 13 | 13,804,071 | 279.7 s | 80,282 |

- In the compiled run the failures are in seven scenarios: the three issue instances
  (`bicritical_cycle_issue502_instance` 87/127, `jacobian_3torsion_genus3_issue574_instance`
  566/732, `rc_circuit_badly_scaled_issue414_instance` 48/144) and the four Fano quintic
  scenarios (8–22 paths each, almost all step size too small). Every other scenario tracks
  all paths.
- Interpreted evaluation rounds differently from compiled: 3.7% more steps, and 142 paths end
  at a different endpoint than the reference: 141 on `jacobian_3torsion_genus3_issue574_instance`
  (singular endpoints) and one on `steiner_conics_real_g2`. The certified count barely moves.
- Scaling invariance: under the random scaling (cᵢ, dⱼ ∈ 10^U(-3, 3)) 193 more paths fail
  than without it (max steps 211 → 357), and 14 fewer distinct solutions are certified.
- `solve` (interpreted, with the endgame): 80,591 successes (80,235 nonsingular, 283
  singular), 25 at infinity, 259 failed, 80,260 distinct certified, 229.5 s. On
  `rc_circuit_badly_scaled_issue414_instance` only 9 of 144 paths succeed, compared with 48
  without the endgame.
- The 27 start sets of the artifact `benchmark_starts` certify complete and distinct (31,078
  solutions, `certify.jl --starts`).
