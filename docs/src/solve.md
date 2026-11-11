# The solve function

The [`solve`](@ref) function is the most convenient way to solve general polynomial systems.
For the mathematical background take a look at our [introduction](https://www.juliahomotopycontinuation.org/guides/introduction/) guide.

```@docs
solve
solver_startsolutions
Solver
```

## Optional path repair

`solve(...; path_repair=true)` runs an aggressive first pass, then retracks failed,
inaccurate, or conflicting regular paths with conservative tracker parameters.
The default is `false`. The option also works when solving for several target
parameter vectors; each target has its own repair pass.

This policy requires a complete batch, so it cannot be combined with
`iterator_only` or `stop_early_cb`.

The function [`paths_to_track`](@ref) allows you to know beforehand how manys the you
need to track:

```@docs
paths_to_track
```
