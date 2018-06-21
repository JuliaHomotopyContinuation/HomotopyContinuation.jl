module Problems

using Compat

import DynamicPolynomials
import MultivariatePolynomials
const MP = MultivariatePolynomials


import ..Input: AbstractInput, TotalDegree, StartTarget, ParameterSystem, MPPolys
import ..Homotopies: AbstractHomotopy, StraightLineHomotopy, ParameterHomotopy
import ..ProjectiveVectors
import ..Systems: AbstractSystem, SPSystem, FPSystem
import ..Systems

using ..Utilities

# STAGE 2 exports
export AbstractProblem,
    AbstractProjectiveProblem,
    AbstractHomogenization,
    NullHomogenization,
    Homogenization,
    Projective,
    homotopy,
    homogenization,
    embed

const DEFAULT_SYSTEM = FPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

abstract type AbstractProblem end


# The HomogenizationStrategy will be useful if we have multi-homogenous stuff
"""
    AbstractHomogenization

Abstract type for describing the homogenization strategy.
"""
abstract type AbstractHomogenization end

"""
    NullHomogenization()

Homogenization was not necessary since it is already homogenous.
"""
struct NullHomogenization <: AbstractHomogenization end

"""
    Homogenization(homvaridx::Int)

Homogenization happend by adding a new variable with index `homvaridx`.
"""
struct Homogenization <: AbstractHomogenization
    homvaridx::Int
end
function Homogenization(homvar::MP.AbstractVariable, variables::Vector{<:MP.AbstractVariable})
    homvaridx = findfirst(var -> var == homvar, variables)
    if homvaridx == 0 || homvaridx === nothing # 0.6 and 0.7
        throw(error("homvar $homvar doesn't occur in the variables!"))
    end
    Homogenization(homvaridx)
end

"""
    Projective(prob::TotalDegreeProblem; options...)
    Projective(prob::StartTargetProblem; options...)
    Projective(start<:Vector{<:Polynomial}, target::Vector{<:MP.Polynomial}; options...)

Construct a `Projective`. This steps constructs a homotopy and homogenizes
the systems if necessary. The options are
* `system=FPSystem`: A constructor to assemble a [`Systems.AbstractSystem`](@ref). The constructor
is called with `system(polys, variables)` where `variables` determines the variable ordering.
* `homotopy=StraightLineHomotopy`: A constructor to construct a [`Homotopies.AbstractHomotopy`](@ref) an `Systems.AbstractSystem`. The constructor
is called with `homotopy(start, target)` where `start` and `target` are systems constructed
with `system`.

    Projective(H::AbstractHomotopy)

Construct a `ProjectiveProblem`. The homotopy `H` needs to be homogenous.
"""
struct Projective{H<:AbstractHomotopy, HS<:AbstractHomogenization} <: AbstractProblem
    homotopy::H
    homogenization::HS
end
function Projective(G::AbstractSystem, F::AbstractSystem,
        homogenization::AbstractHomogenization; homotopy=DEFAULT_HOMOTOPY)
    Projective(homotopy(G, F), homogenization)
end


# We dispatch every input to the two argument case
problem_startsolutions(prob::AbstractInput; homvar=nothing, kwargs...) = problem_startsolutions(prob, homvar; kwargs...)

# TOTALDEGREE
function problem_startsolutions(prob::TotalDegree{Vector{AP}}, _homvar::Nothing; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomial}
    F = prob.system
    variables = MP.variables(F)
    F_ishom = ishomogenous(F)
    if F_ishom
        G = Systems.TotalDegreeSystem(F, variables, variables[1])
        start = totaldegree_solutions(F, NullHomogenization())
        Projective(G, system(F), NullHomogenization(); kwargs...), start
    else
        # We create a new variable to homogenize the system
        homvar = uniquevar(F)
        homogenization = Homogenization(1)
        var_ordering = [homvar; variables]
        proj = Projective(
            Systems.TotalDegreeSystem(F, var_ordering, homvar),
            system(homogenize(F, homvar), var_ordering), homogenization; kwargs...)
        start = totaldegree_solutions(F, homogenization)
        proj, start
    end
end
function problem_startsolutions(prob::TotalDegree{Vector{AP}}, homvar::MP.AbstractVariable; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomialLike}
    @assert ishomogenous(prob.system) "Input system is not homogenous although `homvar=$(homvar)` was passed."
    variables = MP.variables(prob.system)
    homogenization = Homogenization(homvar, variables)
    start = totaldegree_solutions(prob.system, homogenization)

    proj = Projective(
        Systems.TotalDegreeSystem(prob.system, variables, homvar),
        system(homogenize(prob.system, homvar), variables), homogenization; kwargs...)
    proj, start
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    degrees = check_homogenous_degrees(prob.system)
    # system needs to be homogenous
    @assert n + 1 == N "Input system is not a square homogenous system!"

    DynamicPolynomials.@polyvar z[1:N]
    G = Systems.TotalDegreeSystem(degrees, collect(2:N), 1)
    homogenization = NullHomogenization()

    proj = Projective(G, prob.system, homogenization; kwargs...)
    start = totaldegree_solutions(degrees, homogenization)

    proj, start
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Int; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    degrees = check_homogenous_degrees(prob.system)
    # system needs to be homogenous
    @assert n + 1 == N "Input system is not a square homogenous system!"

    homogenization = Homogenization(homvaridx)
    G = Systems.TotalDegreeSystem(degrees, [1:homvaridx-1;homvaridx+1:N], homvaridx)
    proj = Projective(G, prob.system, homogenization; kwargs...)
    start = totaldegree_solutions(degrees, homogenization)

    proj, start
end

# START TARGET
function problem_startsolutions(prob::StartTarget{Vector{AP1}, Vector{AP2}, V}, homvar; system=DEFAULT_SYSTEM, kwargs...) where
    {AP1<:MP.AbstractPolynomialLike, AP2<:MP.AbstractPolynomialLike, V}

    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
    if F_ishom && G_ishom && homvar !== nothing
        Projective(system(G), system(F), Homogenization(homvar, MP.variables(F)); kwargs...),
        prob.startsolutions
    elseif F_ishom && G_ishom && homvar === nothing
        Projective(system(G), system(F), NullHomogenization(); kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        throw(error("One of the input polynomials is homogenous and the other not!"))
    else
        if homvar !== nothing
            throw(error("Input system is not homogenous although `homvar` was passed."))
        end
        homvar = uniquevar(F)
        homogenization = Homogenization(1)
        var_ordering = [homvar; MP.variables(F)]
        Projective(
            system(homogenize(G, homvar), var_ordering),
            system(homogenize(F, homvar), var_ordering), homogenization; kwargs...), prob.startsolutions
    end
end

function problem_startsolutions(prob::ParameterSystem, homvar; system=FPSystem, kwargs...)
    variables = setdiff(allvariables(prob.system), prob.parameters)
    if ishomogenous(prob.system, variables)
        if homvar === nothing
            homogenization = NullHomogenization()
            var_ordering = [variables; prob.parameters]
            var_indices = 1:length(variables)
            param_indices = (length(variables)+1):(length(variables)+length(prob.parameters))
            f = system(prob.system, var_ordering)
        else
            homogenization = Homogenization(homvar, variables)
            var_ordering = [variables; prob.parameters]
            var_indices = 1:length(variables)
            param_indices = (length(variables)+1):(length(variables)+length(prob.parameters))
            f = system(prob.system, var_ordering)
        end
    elseif homvar === nothing
        homogenization = Homogenization(1)
        homvar = uniquevar(prob.system)
        var_ordering = [homvar; variables; prob.parameters]
        var_indices = 1:(length(variables)+1)
        param_indices = (length(variables)+2):(length(variables)+length(prob.parameters)+1)
        f = system(homogenize(prob.system, variables, homvar), var_ordering)
    else
        throw(error("Input system is not homogenous although `homvar` was passed."))
    end

    H = ParameterHomotopy(f, collect(var_indices), collect(param_indices), prob.start, prob.target)

    Projective(H, homogenization), prob.startsolutions
end


function check_homogenous_degrees(F::AbstractSystem)
    n, N = size(F)
    if n ≠ N - 1
        throw(error("Input system is not homogenous! It has $n polynomials in $N variables according to `size`."))
    end
    # The number of variables match, but it still cannot be homogenous.
    # We evaluate the system with y:=rand(N) and 2y. If homogenous then the output
    # scales accordingly to the degrees which we can obtain by taking logarithms.
    x = randn(N)
    cache = Systems.cache(F, x)
    y = Systems.evaluate(F, x, cache)
    scale!(x, 2)
    y2 = Systems.evaluate(F, x, cache)

    degrees = map(y2, y) do y2ᵢ, yᵢ
        # y2ᵢ = 2^dᵢ yᵢ
        float_dᵢ = log2(abs(y2ᵢ / yᵢ))
        dᵢ = round(Int, float_dᵢ)
        if abs(dᵢ - float_dᵢ) > 1e-10
            throw(error("Input system is not homogenous by our numerical check."))
        end
        dᵢ
    end
    degrees
end

"""
    homotopy(prob::Projective)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::Projective) = prob.homotopy

"""
    homogenization(prob::Projective)

Get the [`HomogenizationStrategy`](@ref) stored in the problem `prob`.
"""
homogenization(prob::Projective) = prob.homogenization

"""
    embed(prob::Projective, x)

Embed the solution `x` into projective space if necessary.
"""
function embed(prob::Projective{<:AbstractHomotopy, NullHomogenization}, x)
    M, N = size(homotopy(prob))
    length(x) != N && throw(error("The length of the intial solution is $(length(x)) but expected length $N."))
    return ProjectiveVectors.PVector(x, nothing)
end
function embed(prob::Projective{<:AbstractHomotopy, Homogenization}, x)
    M, N = size(homotopy(prob))
    if length(x) == N
        return ProjectiveVectors.PVector(x, homogenization(prob).homvaridx)
    elseif length(x) == N - 1
        return ProjectiveVectors.embed(x, homogenization(prob).homvaridx)
    end
    throw(error("The length of the initial solution is $(length(x)) but expected length $N or $(N-1)."))
end


"""
    totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})

Returns an iterator of the solutions of the total degree startsystem of `F`.
"""
function totaldegree_solutions(F::Vector{AP}, homogenization) where {AP<:MP.AbstractPolynomialLike}
    totaldegree_solutions(MP.maxdegree.(F), homogenization)
end
function totaldegree_solutions(degrees::Vector{Int}, ::NullHomogenization)
    TotalDegreeSolutionIterator(degrees, true) |> collect
end
function totaldegree_solutions(degrees::Vector{Int}, ::Homogenization)
    TotalDegreeSolutionIterator(degrees, false) |> collect
end


"""
    TotalDegreeSolutionIterator(degrees, homogenous::Bool)

Given the `Vector`s `degrees` `TotalDegreeSolutionIterator` enumerates all solutions
of the system
```math
\\begin{align*}
    z_1^{d_1} - 1 &= 0 \\\\
    z_1^{d_2} - 1 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - 1 &= 0 \\\\
\\end{align*}
```
where ``d_i`` is `degrees[i]`. If `homogenous` is `true` then
the solutions of the system will be embedded in projective space by the map
`x → [1.0; x]`.
"""
struct TotalDegreeSolutionIterator{Iter}
    degrees::Vector{Int}
    homogenous::Bool
    iterator::Iter
end
function TotalDegreeSolutionIterator(degrees::Vector{Int}, homogenous::Bool)
    iterator = Base.Iterators.product(map(d -> 0:d-1, degrees)...)
    TotalDegreeSolutionIterator(degrees, homogenous, iterator)
end

Base.start(iter::TotalDegreeSolutionIterator) = start(iter.iterator)
function Base.next(iter::TotalDegreeSolutionIterator, state)
    indices, nextstate = next(iter.iterator, state)

    value = Vector{Complex{Float64}}()
    if iter.homogenous
        push!(value, complex(1.0, 0.0))
    end
    for (i, k) in zip(1:length(indices), indices)
        push!(value, cis(2π * k / iter.degrees[i]))
    end
    value, nextstate
end
Base.done(iter::TotalDegreeSolutionIterator, state) = done(iter.iterator, state)
Base.length(iter::TotalDegreeSolutionIterator) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeSolutionIterator}) = Vector{Complex{Float64}}

end
