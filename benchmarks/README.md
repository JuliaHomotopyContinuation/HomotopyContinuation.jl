# Benchmarks

This file contains different scripts of polynomial systems and problems we encountered
during the development of this package.
These problems take too much time to include in our normal test suite, but nonetheless
these should be checked manually after major changes in the package.

- `systems/`: a solver-independent collection of benchmark polynomial systems, one TOML
  file per system with provenance and solution count (see `systems/README.md`).
- `harness/`: the path tracking harness, which tracks fixed scenarios built from those
  systems and logs robustness and cost per commit to `log/` (see `harness/README.md`).
