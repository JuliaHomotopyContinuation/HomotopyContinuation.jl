# Homotopies

## Pre-defined homotopies
```@docs
StraightLineHomotopy
FixedPointHomotopy
```
## Interface for custom homotopies

### Abstract types
```@docs
Homotopies.AbstractHomotopy
Homotopies.AbstractHomotopyCache
```

### Mandatory
```@docs
Homotopies.evaluate!
Homotopies.jacobian!
Homotopies.dt!
Base.size(::Homotopies.AbstractHomotopy)
```
### Optional
```@docs
Homotopies.cache
Homotopies.evaluate_and_jacobian!
Homotopies.evaluate_and_jacobian
Homotopies.jacobian_and_dt!
Homotopies.evaluate
Homotopies.jacobian
Homotopies.dt
Homotopies.precondition!
Homotopies.update!
```
