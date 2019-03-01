# Newton's method

Sometimes it is necessary to refine obtained solutions.
For this we provide an interface to Newton's method.

```@docs
newton
NewtonResult
```

For high performance applications we also provide an in-place version of Newton's method
which avoids any temporary allocations.
```@docs
newton!
NewtonCache
```
