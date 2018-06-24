using Compat

import MultivariatePolynomials
const MP = MultivariatePolynomials

import .Input
import .Problems
import .Systems
import .Homotopies
import .Solving

using .Utilities

export solve


"""
    solve(F; options...)

Solve the system `F` using a total degree homotopy. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- [`Systems.AbstractSystem`](@ref) (the system has to represent a **homogenous** polynomial system.)

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

To solve ``F`` by a custom `Systems.AbstractSystem` you can do
```julia
@polyvar x y z
# The system `F` has to be homgoenous system
F = Systems.SPSystem([x^2+y^2+z^2, 2x+3y-z]) # Systems.SPSystem <: Systems.AbstractSystem
# To solve the original affine system we have to tell that the homogenization variable has index 3
solve(F, homvar=3)
```
or equivalently (in this case) by
```julia
solve([x^2+y^2+z^2, 2x+3y-z], system=Systems.SPSystem)
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
    solve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial}, parametervariables, startparameters, targetparameters, startsolutions)

Construct a parameter homotopy ``H(x,t) := F(x; tp₁+(1-t)p₀)`` where `p₁` is `startparameters`, `p₀` is `targetparameters`
and the `parametervariables` are the variables of `F` which should be considerd parameters.


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
p₁ = [1, 0]
p₀ = [2, 4]
startsolutions = [[1, 1]]
solve(F, a, p₁, p₀, startsolutions)
```

# Abstract Homotopy

    solve(H::Homotopies.AbstractHomotopy, start_solutions; options...)

Solve the homotopy `H` by tracking the each solution of
``H(⋅, t)`` (as provided by `start_solutions`) from ``t=1`` to ``t=0``.
Note that `H` has to be a homotopy between *homogenous* polynomial systems.
If it should be considered as an affine system indicate which is the index
of the homogenization variable, e.g. `solve(H, startsolutions, homvar=3)`
if the third variable is the homogenization variable.


# Options
General options:

* `system::Systems.AbstractSystem`: A constructor to assemble a [`Systems.AbstractSystem`](@ref). The default is [`Systems.FPSystem`](@ref). This constructor is only applied to the input of `solve`. The constructor is called with `system(polynomials, variables)` where `polynomials` is a vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the variable ordering.
* `homotopy::Systems.AbstractHomotopy`: A constructor to construct a [`Homotopies.AbstractHomotopy`](@ref). The default is [`StraightLineHomotopy`](@ref). The constructor is called with `homotopy(start, target)` where `start` and `target` are homogenous [`Systems.AbstractSystem`](@ref)s.
* `seed::Int`: The random seed used during the computations.
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the *homogenous* system `F` as an affine system which was homogenized by `homvar`. If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`) `homvar` is the actual variable used in the system `F`.
* `endgame_start=0.1`: The value of `t` for which the endgame is started.
* `report_progress=true`: Whether a progress bar should be printed to `STDOUT`.

Pathtracking specific:
* `corrector::Correctors.AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`Correctors.Newton`](@ref).
* `corrector_maxiters=2`: The maximal number of correction steps in a single step.
* `predictor::Predictors.AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Predictors.RK4`](@ref).
* `refinement_maxiters=corrector_maxiters`: The maximal number of correction steps used to refine the final value.
* `refinement_tol=1e-11`: The precision used to refine the final value.
* `tol=1e-7`: The precision used to track a value.
* `initial_steplength=0.1`: The initial step size for the predictor.
* `steplength_increase_factor=2.0`: The factor with which the step size is increased after `steplength_consecutive_successes_necessary` consecutive successes.
* `steplength_decrease_factor=inv(increase_factor)`: The factor with which the step size is decreased after a step failed.
* `steplength_consecutive_successes_necessary=5`: The numer of consecutive successes necessary until the step size is increased by `steplength_increase_factor`.
* `maximal_steplength=max(0.1, initial_steplength)`: The maximal step length.
* `minimal_steplength=1e-14`: The minimal step size. If the size of step is below this the path is considered failed.

Endgame specific options
* `cauchy_loop_closed_tolerance=1e-3`: The tolerance for which is used to determine whether a loop is closed. The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.
* `cauchy_samples_per_loop=6`: The number of samples used to predict an endpoint. A higher number of samples should result in a better approximation. Note that the error should be roughly ``t^n`` where ``t`` is the current time of the loop and ``n`` is `cauchy_samples_per_loop`.
* `egtol=1e-10`: This is the tolerance necessary to declare the endgame converged.
* `maxnorm=1e5`: If our original problem is affine we declare a path at infinity if the infinity norm with respect to the standard patch is larger than `maxnorm`.
* `maxwindingnumber=15`: The maximal windingnumber we try to find using Cauchys integral formula.
* `max_extrapolation_samples=4`: During the endgame a Richardson extrapolation is used to improve the accuracy of certain approximations. This is the maximal number of samples used for this.
* `minradius=1e-15`: A path is declared false if the endgame didn't finished until then.
* `sampling_factor=0.5`: During the endgame we approach ``0`` by the geometric series ``h^kR₀`` where ``h`` is `sampling_factor` and `R₀` the endgame start provided in `runendgame`.

"""
function solve end

# External
function solve(F::Vector{<:MP.AbstractPolynomial}; seed=randseed(), homvar=nothing, kwargs...)
    srand(seed)
    F = filter(f -> !iszero(f), F)
    checkfinite_dimensional(F, homvar)
    solve(Input.TotalDegree(F), seed; homvar=homvar, kwargs...)
end

function solve(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions; homvar=nothing, seed=randseed(), kwargs...)
    srand(seed)
    @assert length(G) == length(F)
    checkfinite_dimensional(F, homvar)
    solve(Input.StartTarget(G, F, promote_startsolutions(startsolutions)), seed; homvar=homvar, kwargs...)
end

function solve(F::Systems.AbstractSystem; seed=randseed(), kwargs...)
    srand(seed)
	solve(Input.TotalDegree(F), seed; kwargs...)
end

function checkfinite_dimensional(F::Vector{<:MP.AbstractPolynomial}, homvar)
    N = homvar === nothing ? MP.nvariables(F) : MP.nvariables(F) - 1
    n = length(F)
    # square system and each polynomial is non-zero
    if n ≥ N ||
       n == N - 1 && ishomogenous(F)
        return
    end
    throw(AssertionError("The input system will not result in a finite number of solutions."))
end

function solve(F::Vector{<:MP.AbstractPolynomial},
    p::Vector{<:MP.AbstractVariable},
    a_1::Vector{<:Number},
    a_2::Vector{<:Number},
    startsolutions; seed=randseed(), homotopy=nothing, kwargs...)
    srand(seed)

    @assert length(p) == length(a_1) "Number of parameters must match"
    @assert length(a_1) == length(a_2) "Start and target parameters must have the same length"

    STP =
    solve(Input.ParameterSystem(F, p, a_1, a_2, promote_startsolutions(startsolutions)), seed; kwargs...)
end

function solve(H::Homotopies.AbstractHomotopy, startsolutions; seed=randseed(), kwargs...)
    srand(seed)
	solve(Input.Homotopy(H, promote_startsolutions(startsolutions)), seed; kwargs...)
end

promote_startsolutions(xs::Vector{Vector{Complex128}}) = xs
function promote_startsolutions(xs::Vector{<:AbstractVector{<:Number}})
    PT = promote_type(typeof(xs[1][1]), Complex{Float64})
    map(s -> convert.(PT, s), xs)
end

randseed() = rand(1_000:1_000_000)

# Internal
function solve(input::Input.AbstractInput, seed;
	homvar::Union{Nothing, Int, MP.AbstractVariable}=nothing,
    system=Systems.FPSystem,
    homotopy=Homotopies.StraightLineHomotopy,
    kwargs...)

    p, startsolutions = Problems.problem_startsolutions(input, homvar=homvar, system=system, homotopy=homotopy)
    solve(p, startsolutions, seed; kwargs...)
end

function solve(prob::Problems.AbstractProblem, start_solutions, seed; threading=true, kwargs...)
    solver = Solving.Solver(prob, start_solutions, 1.0, 0.0, seed; kwargs...)
    if threading && Threads.nthreads() > 1
        solvers = append!([solver], [deepcopy(solver) for _=2:Threads.nthreads()])
        Solving.solve(solvers, start_solutions)
    else
        Solving.solve(solver, start_solutions)
    end
end
