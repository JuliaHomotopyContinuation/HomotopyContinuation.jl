# Monodromy Solve

Next to [`solve`](@ref), HomotopyContinuation.jl provides the function `monodromy_solve`. Instead of taking a two systems `f` and `g` and tracking an array of start solutions from `f` to `g`, `monodromy_solve` takes as input a single system with parameters `p` and together with a start solution `s`. Then by tracking `s` around loops in the parameters `p`, `monodromy_solve` duplicates the solution until a stopping criterion is reached.

```@docs
monodromy_solve
```

## Strategies
```@docs
Triangle
Petal
```

## GroupActions

If there is a group acting on the solution set of the polynomial system this can provided with the `group_action` keyword for single group actions or with the `group_actions` keyword for compositions
of group actions. These will be internally transformed into `GroupActions`.

```@docs
GroupActions
```

To help with the more common group actions we provide some helper functions:

```@docs
SymmetricGroup
```


## Helper functions
```@docs
parameters
```
