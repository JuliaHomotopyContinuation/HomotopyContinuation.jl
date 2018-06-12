# Polynomial systems

## Pre-defined systems
```@docs
FPSystem
SPSystem
```

## Interface for custom systems

### Abstract types
```@docs
Systems.AbstractSystem
Systems.AbstractSystemCache
```

### Mandatory
```@docs
Systems.evaluate!
Systems.evaluate
Systems.jacobian!
Systems.jacobian
Base.size(::Systems.AbstractSystem)
```
### Optional
```@docs
Systems.cache
Systems.evaluate_and_jacobian!
Systems.evaluate_and_jacobian
```
