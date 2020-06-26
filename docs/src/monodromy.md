# Solving parametrized systems with monodromy

Next to [`solve`](@ref), HomotopyContinuation.jl provides the function
[`monodromy_solve`](@ref) which uses the [monodromy method](https://www.juliahomotopycontinuation.org/guides/monodromy/#monodromy)
to solve a parameterized system of polynomials.
Often [`monodromy_solve`](@ref) allows to still compute all isolated solutions of system where
the number of paths tracked in  `solve`](@ref) is already infeasible.
Make sure to check out our
[monodromy guide](https://www.juliahomotopycontinuation.org/guides/monodromy/)
for a more in depth introduction into this method.

```@docs
monodromy_solve
MonodromyOptions
find_start_pair
```

## Monodromy Result

A call to [`monodromy_solve`](@ref) returns a [`MonodromyResult`](@ref):
```@docs
MonodromyResult
solutions(::MonodromyResult)
nsolutions(result::MonodromyResult)
parameters(::MonodromyResult)
results(result::MonodromyResult)
nresults(result::MonodromyResult)
is_success(result::MonodromyResult)
is_heuristic_stop(result::MonodromyResult)
seed(result::MonodromyResult)
```

## Group actions

If there is a group acting on the solution set of the polynomial system this can provided with the `group_action` keyword for single group actions or with the `group_actions` keyword for compositions
of group actions. These will be internally transformed into `GroupActions`.

```@docs
GroupActions
```

To help with the more common group actions we provide some helper functions:

```@docs
SymmetricGroup
```
