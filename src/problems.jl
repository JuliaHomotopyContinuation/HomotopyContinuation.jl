module Problems

import DynamicPolynomials
import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials
import Random

import ..Input
import ..Input: AbstractInput, TotalDegree, StartTarget, ParameterSystem, MPPolys
import ..Homotopies: AbstractHomotopy, StraightLineHomotopy, ParameterHomotopy
import ..ProjectiveVectors
import ..ProjectiveVectors: PVector
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
    Projective(H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `ProjectiveProblem`. The homotopy `H` needs to be homogenous.
"""
struct Projective{H<:AbstractHomotopy, HS<:AbstractHomogenization} <: AbstractProblem
    homotopy::H
    homogenization::HS
    seed::Int
end
function Projective(G::AbstractSystem, F::AbstractSystem,
        homogenization::AbstractHomogenization, seed::Int; homotopy=DEFAULT_HOMOTOPY)
    Projective(homotopy(G, F), homogenization, seed)
end

const supported_kwargs = [:seed, :homvar, :homotopy, :system]

Base.broadcastable(P::AbstractProblem) = Ref(P)
# We dispatch every input to the two argument case
"""
    problem_startsolutions(F; options...)
    problem_startsolutions(G, F, startsolutions; options...)
    problem_startsolutions(H::AbstractHomotopy, startsolutions; options...)
    problem_startsolutions(prob::TotalDegreeProblem; options...)
    problem_startsolutions(prob::StartTargetProblem; options...)

    Construct the problem and if necessary the startsolutions. This steps
    constructs a homotopy and homogenizes the systems if necessary.

    The `options` are
    * `seed::Int`: Random seed used in the construction.
    * `system=FPSystem`: A constructor to assemble a [`Systems.AbstractSystem`](@ref). The constructor
    is called with `system(polys, variables)` where `variables` determines the variable ordering.
    * `homotopy=StraightLineHomotopy`: A constructor to construct a [`Homotopies.AbstractHomotopy`](@ref) an `Systems.AbstractSystem`. The constructor
    is called with `homotopy(start, target)` where `start` and `target` are systems constructed
    with `system`.
"""
function problem_startsolutions end

function problem_startsolutions(args...; seed=randseed(), kwargs...)
    Random.seed!(seed)
    _problem_startsolutions(args..., seed; kwargs...)
end


#########
# Input
#########

function _problem_startsolutions(F::Vector{<:MP.AbstractPolynomial}, seed::Int; homvar=nothing, kwargs...)
    F = filter(f -> !iszero(f), F)
    check_zero_dimensional(F, homvar)
    # square system and each polynomial is non-zero
    if length(F) == MP.nvariables(F) && ishomogenous(F)
        throw(AssertionError("The input system is a square homogenous system. This will result in an at least 1 dimensional solution space."))
    end
    _problem_startsolutions(Input.TotalDegree(F), seed; homvar=homvar, kwargs...)
end

function _problem_startsolutions(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions, seed; homvar=nothing, kwargs...)
    @assert length(G) == length(F) "Start and target system don't have the same length"
    check_zero_dimensional(F, homvar)
    _problem_startsolutions(Input.StartTarget(G, F, promote_startsolutions(startsolutions)), seed; homvar=homvar, kwargs...)
end

function _problem_startsolutions(F::Systems.AbstractSystem, seed; kwargs...)
	_problem_startsolutions(Input.TotalDegree(F), seed; kwargs...)
end

function _problem_startsolutions(F::Vector{<:MP.AbstractPolynomial},
    p::Vector{<:MP.AbstractVariable},
    a_1::Vector{<:Number},
    a_2::Vector{<:Number},
    startsolutions, seed::Int; homotopy=nothing, kwargs...)

    @assert length(p) == length(a_1) "Number of parameters must match"
    @assert length(a_1) == length(a_2) "Start and target parameters must have the same length"

    _problem_startsolutions(Input.ParameterSystem(F, p, a_1, a_2, promote_startsolutions(startsolutions)), seed; kwargs...)
end

function _problem_startsolutions(H::AbstractHomotopy, startsolutions, seed::Int; kwargs...)
    _problem_startsolutions(Input.Homotopy(H, promote_startsolutions(startsolutions)), seed; kwargs...)
end

function _problem_startsolutions(input::AbstractInput, seed::Int;
	homvar::Union{Nothing, Int, MP.AbstractVariable}=nothing, kwargs...)
    _problem_startsolutions(input, homvar, seed; kwargs...)
end

function _problem_startsolutions(input::Input.Homotopy, homvar::Int, seed; kwargs...)
    check_homogenous_degrees(Systems.FixedHomotopy(input.H, rand()))
    Projective(input.H, Homogenization(homvar), seed), input.startsolutions
end

function _problem_startsolutions(input::Input.Homotopy, homvar::Nothing, seed; kwargs...)
    check_homogenous_degrees(Systems.FixedHomotopy(input.H, rand()))
    Projective(input.H, NullHomogenization(), seed), input.startsolutions
end


promote_startsolutions(iter) = promote_startsolutions(collect(iter))
promote_startsolutions(xs::Vector{Vector{ComplexF64}}) = xs
function promote_startsolutions(xs::Vector{<:AbstractVector{<:Number}})
    PT = promote_type(typeof(xs[1][1]), Complex{Float64})
    map(s -> convert.(PT, s), xs)
end


const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""


##############
# TOTALDEGREE
##############

function _problem_startsolutions(prob::TotalDegree{Vector{AP}}, _homvar::Nothing, seed::Int; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomial}
    F, variables, homogenization = homogenize_if_necessary(prob.system)
    if homogenization == NullHomogenization()
        G = Systems.TotalDegreeSystem(F, variables, variables[1])
        start = totaldegree_solutions(F, NullHomogenization())
        Projective(G, system(F), NullHomogenization(), seed; kwargs...), start
    elseif homogenization isa Homogenization
        proj = Projective(
            Systems.TotalDegreeSystem(F, variables, variables[homogenization.homvaridx]),
            system(F, variables), homogenization, seed; kwargs...)
        proj, totaldegree_solutions(F, homogenization)
    end
end

function _problem_startsolutions(prob::TotalDegree{Vector{AP}},
    homvar::MP.AbstractVariable, seed; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomialLike}
    @assert ishomogenous(prob.system) "Input system is not homogenous although `homvar=$(homvar)` was passed."
    F, variables, homogenization = homogenize_if_necessary(prob.system; homvar=homvar)

    start = totaldegree_solutions(prob.system, homogenization)
    proj = Projective(
        Systems.TotalDegreeSystem(prob.system, variables, homvar),
        system(F, variables), homogenization, seed; kwargs...)
    proj, start
end

function _problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    degrees = check_homogenous_degrees(prob.system)
    # system needs to be homogenous
    if n + 1 > N
        throw(AssertionError(overdetermined_error_msg))
    end
    @assert n + 1 == N "Input system is not a square homogenous system!"

    DynamicPolynomials.@polyvar z[1:N]
    G = Systems.TotalDegreeSystem(degrees, collect(2:N), 1)
    homogenization = NullHomogenization()

    proj = Projective(G, prob.system, homogenization, seed; kwargs...)
    start = totaldegree_solutions(degrees, homogenization)

    proj, start
end

function _problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Int, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    # system needs to be homogenous
    degrees = check_homogenous_degrees(prob.system)
    if n + 1 > N
        throw(AssertionError(overdetermined_error_msg))
    end
    @assert n + 1 == N "Input system is not a square homogenous system!"

    homogenization = Homogenization(homvaridx)
    G = Systems.TotalDegreeSystem(degrees, [1:homvaridx-1;homvaridx+1:N], homvaridx)
    proj = Projective(G, prob.system, homogenization, seed; kwargs...)
    start = totaldegree_solutions(degrees, homogenization)

    proj, start
end


###############
# START TARGET
###############

function _problem_startsolutions(prob::StartTarget{Vector{AP1}, Vector{AP2}, V}, homvar, seed; system=DEFAULT_SYSTEM, kwargs...) where
    {AP1<:MP.AbstractPolynomialLike, AP2<:MP.AbstractPolynomialLike, V}

    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
    if F_ishom && G_ishom && homvar !== nothing
        Projective(system(G), system(F), Homogenization(homvar, MP.variables(F)), seed; kwargs...),
        prob.startsolutions
    elseif F_ishom && G_ishom && homvar === nothing
        Projective(system(G), system(F), NullHomogenization(), seed; kwargs...), prob.startsolutions
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
            system(homogenize(F, homvar), var_ordering), homogenization, seed; kwargs...), prob.startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function _problem_startsolutions(prob::ParameterSystem, homvar, seed; system=FPSystem, kwargs...)
    F, variables, homogenization = homogenize_if_necessary(prob.system, homvar=homvar, parameters=prob.parameters)

    H = ParameterHomotopy(F, variables, prob.parameters, prob.p₁, prob.p₀)

    Projective(H, homogenization, seed), prob.startsolutions
end

##########
# HELPERS
##########

"""
    homogenize_if_necessary(F::Vector{<:MP.AbstractPolynomialLike})

Homogenizes the system `F` if necessary and returns the (new) system `F` its variables
and a subtype of [`AbstractHomogenization`] indicating whether it was homegenized.
If it was homogenized and no then the new variable is the **first one**.
"""
function homogenize_if_necessary(F::Vector{<:MP.AbstractPolynomialLike}; homvar=nothing, parameters=nothing)
    variables = MP.variables(F)
    if parameters !== nothing
        variables = setdiff(variables, parameters)
    end

    n, N = length(F), length(variables)
    if ishomogenous(F; parameters=parameters)
        # N = n+1 is the only valid size configuration
        if n + 1 > N
            throw(AssertionError(overdetermined_error_msg))
        end
        if homvar === nothing
            F, variables, NullHomogenization()
        else
            F, variables, Homogenization(homvar, variables)
        end
    else
        # N = n is the only valid size configuration
        if n > N
            throw(AssertionError(overdetermined_error_msg))
        end
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end
        # We create a new variable to homogenize the system
        homvar = uniquevar(F)
        push!(variables, homvar)
        sort!(variables, rev=true)

        homogenize(F, homvar; parameters=parameters), variables, Homogenization(homvar, variables)
    end
end

function check_homogenous_degrees(F::AbstractSystem)
    n, N = size(F)
    if n < N - 1
        throw(error("Input system is not homogenous! It has $n polynomials in $N variables according to `size`."))
    end
    # The number of variables match, but it still cannot be homogenous.
    # We evaluate the system with y:=rand(N) and 2y. If homogenous then the output
    # scales accordingly to the degrees which we can obtain by taking logarithms.
    x = rand(ComplexF64, N)
    cache = Systems.cache(F, x)
    y = Systems.evaluate(F, x, cache)
    LinearAlgebra.rmul!(x, 2)
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
    length(x) != N && throw(error("The length of the initial solution is $(length(x)) but expected length $N."))
    return PVector(x, nothing)
end
function embed(prob::Projective{<:AbstractHomotopy, Homogenization}, x)
    M, N = size(homotopy(prob))
    if length(x) == N
        return PVector(x, homogenization(prob).homvaridx)
    elseif length(x) == N - 1
        return ProjectiveVectors.embed(x, homogenization(prob).homvaridx)
    end
    throw(error("The length of the initial solution is $(length(x)) but expected length $N or $(N-1)."))
end
embed(prob::Projective{<:AbstractHomotopy, Homogenization}, x::PVector) = x
embed(prob::Projective{<:AbstractHomotopy, NullHomogenization}, x::PVector) = x


"""
    totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})

Returns an iterator of the solutions of the total degree startsystem of `F`.
"""
function totaldegree_solutions(F::Vector{AP}, homogenization) where {AP<:MP.AbstractPolynomialLike}
    totaldegree_solutions(MP.maxdegree.(F), homogenization)
end
function totaldegree_solutions(degrees::Vector{Int}, ::NullHomogenization)
    TotalDegreeSolutionIterator(degrees, true)
end
function totaldegree_solutions(degrees::Vector{Int}, ::Homogenization)
    TotalDegreeSolutionIterator(degrees, false)
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

function Base.iterate(iter::TotalDegreeSolutionIterator)
    indices, state = iterate(iter.iterator)
    _value(iter, indices), state
end
function Base.iterate(iter::TotalDegreeSolutionIterator, state)
    it = iterate(iter.iterator, state)
    it === nothing && return nothing
    _value(iter, it[1]), it[2]
end

function _value(iter::TotalDegreeSolutionIterator, indices)
    value = Vector{Complex{Float64}}()
    if iter.homogenous
        push!(value, complex(1.0, 0.0))
    end
    for (i, k) in enumerate(indices)
        push!(value, cis(2π * k / iter.degrees[i]))
    end
    value
end

Base.length(iter::TotalDegreeSolutionIterator) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeSolutionIterator}) = Vector{Complex{Float64}}

end
