export solve


"""
    solve(F; options...)

Solve the system `F` using a total degree homotopy. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- [`AbstractSystem`](@ref) (the system has to represent a **homogenous** polynomial system.)

### Example
Assume we want to solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)``.
```julia
@polyvar x y
solve([x^2+y^2+1, 2x+3y-1])
```
If you polynomial system is already homogenous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This would result in the same result as `solve([x^2+y^2+1, 2x+3y-1])`.

To solve ``F`` by a custom `AbstractSystem` you can do
```julia
@polyvar x y z
# The system `F` has to be homgoenous system
F = SPSystem([x^2+y^2+z^2, 2x+3y-z]) # SPSystem <: AbstractSystem
# To solve the original affine system we have to tell that the homogenization variable has index 3
solve(F, homvar=3)
```
or equivalently (in this case) by
```julia
solve([x^2+y^2+z^2, 2x+3y-z], system=SPSystem)
```

# Start Target Homotopy

    solve(G, F, start_solutions; options...)

Solve the system `F` by tracking the each provided solution of
`G` (as provided by `start_solutions`).

### Example
```julia
@polyvar x y
G = [x^2-1,y-1]
F = [x^2+y^2+z^2, 2x+3y-z]
solve(G, F, [[1, 1], [-1, 1]])
```

# Parameter Homotopy
    solve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},
        startsolutions; parameters::Vector{<:MP.AbstractVariable}, p₁, p₀, γ₁=nothing, γ₀=nothing)

Solve the parameter homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀))
```,
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input `parameters` specifies the parameter variables of `F`
which should be considered as parameters.
Neccessarily, ``length(parameters) == length(p₁) == length(p₀)``.

    solve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},
            startsolutions; parameters::Vector{<:MP.AbstractVariable},
            startparameters, targetparameters,
            startgamma=randn(ComplexF64), targetgamma=randn(ComplexF64))

This is a non-unicode variant where `γ₁=start_parameters`, `γ₀=target_parameters`,
    `γ₁=start_gamma`, γ₀=`target_gamma`.

## Example
We want to solve a parameter homotopy ``H(x,t) := F(x; t[1, 0]+(1-t)[2, 4])`` where
```math
F(x; a) := (x₁^2-a₁, x₁x₂-a₁+a₂)
```
and let's say we are only intersted in tracking of ``[1,1]``.
This can be accomplished as follows
```julia
@polyvar x[1:2] a[1:2]
F = [x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]]
startsolutions = [[1, 1]]
solve(F, startsolutions, parameters=a, p₁=p₁, p₀=p₀)
# If you don't like unicode this is also possible
solve(F, startsolutions, parameters=a, startparameters=p₁, targetparameters=p₀)
```

# Abstract Homotopy

    solve(H::AbstractHomotopy, start_solutions; options...)

Solve the homotopy `H` by tracking the each solution of
``H(⋅, t)`` (as provided by `start_solutions`) from ``t=1`` to ``t=0``.
Note that `H` has to be a homotopy between *homogenous* polynomial systems.
If it should be considered as an affine system indicate which is the index
of the homogenization variable, e.g. `solve(H, startsolutions, homvar=3)`
if the third variable is the homogenization variable.


# Options
General options:

* `system::AbstractSystem`: A constructor to assemble a [`AbstractSystem`](@ref). The default is [`SPSystem`](@ref). This constructor is only applied to the input of `solve`. The constructor is called with `system(polynomials, variables)` where `polynomials` is a vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the variable ordering. If you experience significant compilation times, consider to change system to `FPSystem`.
* `homotopy::AbstractHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref). The default is [`StraightLineHomotopy`](@ref). The constructor is called with `homotopy(start, target)` where `start` and `target` are homogenous [`AbstractSystem`](@ref)s.
* `seed::Int`: The random seed used during the computations.
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the *homogenous* system `F` as an affine system which was homogenized by `homvar`. If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`) `homvar` is the actual variable used in the system `F`.
* `endgame_start=0.1`: The value of `t` for which the endgame is started.
* `report_progress=true`: Whether a progress bar should be printed to `STDOUT`.
* `threading=true`: Enable or disable multi-threading.

Pathtracking specific:
* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `max_corrector_iters=2`: The maximal number of correction steps in a single step.
* `accuracy=1e-7`: The accuracy used to track a value.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref).
* `max_refinement_iters=max_corrector_iters`: The maximal number of correction steps used to refine the final value.
* `refinement_accuracy=1e-8`: The precision used to refine the final value.
* `initial_step_size=0.1`: The initial step size for the predictor.
* `min_step_size=1e-14`: The minimal step size. If the size of step is below this the path is considered failed.
* `max_steps=1000`: The maximal number of steps per path.

Endgame specific options
* `cauchy_loop_closed_tolerance=1e-3`: The tolerance for which is used to determine whether a loop is closed. The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.
* `cauchy_samples_per_loop=6`: The number of samples used to predict an endpoint. A higher number of samples should result in a better approximation. Note that the error should be roughly ``t^n`` where ``t`` is the current time of the loop and ``n`` is `cauchy_samples_per_loop`.
* `egtol=1e-10`: This is the tolerance necessary to declare the endgame converged.
* `maxnorm=1e5`: If our original problem is affine we declare a path at infinity if the infinity norm with respect to the standard patch is larger than `maxnorm`.
* `maxwindingnumber=15`: The maximal windingnumber we try to find using Cauchys integral formula.
* `max_extrapolation_samples=4`: During the endgame a Richardson extrapolation is used to improve the accuracy of certain approximations. This is the maximal number of samples used for this.
* `minradius=1e-15`: A path is declared false if the endgame didn't finished until then.
* `sampling_factor=0.5`: During the endgame we approach ``0`` by the geometric series ``h^kR₀`` where ``h`` is `sampling_factor` and `R₀` the endgame start provided in `runendgame`.
* `maxiters_per_step=100`: The maximal number of steps bewtween two samples.
"""
function solve end

function solve(args...; threading=true, kwargs...)
    solver, startsolutions = solver_startsolutions(args...; kwargs...)
    solve(solver, startsolutions, threading=threading)
end

# Internal
function solve(solver::Solver, start_solutions; threading=true)
    if threading && Threads.nthreads() > 1
        solvers = append!([solver], [deepcopy(solver) for _=2:Threads.nthreads()])
        solve(solvers, start_solutions, threading=true)
    else
        internal_solve(solver, start_solutions)
    end
end
function solve(solvers::AbstractVector{<:Solver}, start_solutions; threading=true)
    internal_solve(threading ? solvers : solvers[1], start_solutions)
end
