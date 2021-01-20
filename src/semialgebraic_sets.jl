export SemialgebraicSetsHCSolver

function ModelKit.System(V::SemialgebraicSets.AbstractAlgebraicSet; kwargs...)
    ModelKit.System(SemialgebraicSets.equalities(V); kwargs...)
end

"""
    SemialgebraicSetsHCSolver(; options...)

Construct a `SemialgebraicSets.AbstractAlgebraicSolver` to be used in `SemialgebraicSets`.
`options` are all valid options for [`solve`](@ref).

## Example

```
julia> using HomotopyContinuation, SemialgebraicSets;

julia> solver = SemialgebraicSetsHCSolver(; compile = false)
SemialgebraicSetsHCSolver(; compile = false)

julia> @polyvar x y
(x, y)

julia> V = @set x^2 == 1 && y^2 == 2 solver
Algebraic Set defined by 2 equalities
 x^2 - 1.0 = 0
 y^2 - 2.0 = 0

julia> collect(V)
4-element Array{Array{Float64,1},1}:
 [1.0, 1.414213562373095]
 [1.0, -1.414213562373095]
 [-1.0, 1.414213562373095]
 [-1.0, -1.414213562373095]
```
"""
struct SemialgebraicSetsHCSolver <: SemialgebraicSets.AbstractAlgebraicSolver
    options::Any
end
SemialgebraicSetsHCSolver(; options...) = SemialgebraicSetsHCSolver(options)
function Base.show(io::IO, solver::SemialgebraicSetsHCSolver)
    print(io, "SemialgebraicSetsHCSolver(; ")
    join(
        io,
        ["$(k) = $(v isa Symbol ? Expr(:quote, v) : v)" for (k, v) in solver.options],
        ", ",
    )
    print(io, ")")
end

function SemialgebraicSets.solvealgebraicequations(
    V::SemialgebraicSets.AbstractAlgebraicSet,
    solver::SemialgebraicSetsHCSolver,
)
    F = System(V)
    m, n = size(F)
    # Check that we can have isolated solutions
    m ≥ n || return nothing
    # Only return real, non-singular solutions
    return real_solutions(solve(F; solver.options...); only_nonsingular = true)
end
