export problem_startsolutions

const supported_keywords = [[:seed, :homvar, :homotopy, :system]; Input.supported_keywords]
const DEFAULT_SYSTEM = FPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

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
    supported, rest = Utilities.splitkwargs(kwargs, Input.supported_keywords)
    input = Input.input(args...; supported...)
    problem_startsolutions(input, seed; rest...)
end

function problem_startsolutions(input::AbstractInput; seed=randseed(), kwargs...)
    problem_startsolutions(input, seed; kwargs...)
end
function problem_startsolutions(input::AbstractInput, seed::Int;
	homvar::Union{Nothing, Int, MP.AbstractVariable}=nothing, kwargs...)
    problem_startsolutions(input, homvar, seed; kwargs...)
end

function problem_startsolutions(input::Input.Homotopy, homvar, seed; kwargs...)
    Projective(input.H, homogenization(homvar), seed), input.startsolutions
end


# promote_startsolutions(iter) = promote_startsolutions(collect(iter))
# promote_startsolutions(xs::Vector{Vector{ComplexF64}}) = xs
# function promote_startsolutions(xs::Vector{<:AbstractVector{<:Number}})
#     PT = promote_type(typeof(xs[1][1]), Complex{Float64})
#     map(s -> convert.(PT, s), xs)
# end


const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""


##############
# TOTALDEGREE
##############

function problem_startsolutions(prob::TotalDegree{Vector{AP}}, _homvar::Nothing, seed::Int; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomial}
    F, variables, homogenization = homogenize_if_necessary(prob.system)
    if ishomogenized(homogenization)
        proj = Projective(
            Systems.TotalDegreeSystem(prob.degrees, variables, variables[homogenization.homvaridx]),
            system(F, variables), homogenization, seed; kwargs...)
        proj, totaldegree_solutions(prob.degrees, homogenization)
    else
        G = Systems.TotalDegreeSystem(prob.degrees, variables, variables[1])
        start = totaldegree_solutions(prob.degrees, NullHomogenization())
        Projective(G, system(F), NullHomogenization(), seed; kwargs...), start
    end
end

function problem_startsolutions(prob::TotalDegree{Vector{AP}},
    homvar::MP.AbstractVariable, seed; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomialLike}

    if !ishomogenous(prob.system)
        error("Input system is not homogenous although `homvar=$(homvar)` was passed.")
    end
    F, variables, homogenization = homogenize_if_necessary(prob.system; homvar=homvar)

    start = totaldegree_solutions(prob.degrees, homogenization)
    proj = Projective(
        Systems.TotalDegreeSystem(prob.degrees, variables, homvar),
        system(F, variables), homogenization, seed; kwargs...)
    proj, start
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    G = Systems.TotalDegreeSystem(prob.degrees, collect(2:N), 1)

    (Projective(G, prob.system, NullHomogenization(), seed; kwargs...),
     totaldegree_solutions(prob.degrees, NullHomogenization()))
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Int, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)

    homogenization = Homogenization(homvaridx)
    G = Systems.TotalDegreeSystem(prob.degrees, [1:homvaridx-1;homvaridx+1:N], homvaridx)

    (Projective(G, prob.system, homogenization, seed; kwargs...),
     totaldegree_solutions(prob.degrees, homogenization))
end


###############
# START TARGET
###############

function problem_startsolutions(prob::StartTarget{Vector{AP1}, Vector{AP2}}, homvar, seed; system=DEFAULT_SYSTEM, kwargs...) where
    {AP1<:MP.AbstractPolynomialLike, AP2<:MP.AbstractPolynomialLike}

    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
    if F_ishom && G_ishom && homvar !== nothing
        Projective(system(G), system(F), Homogenization(homvar, MP.variables(F)), seed; kwargs...),
        prob.startsolutions
    elseif F_ishom && G_ishom && homvar === nothing
        Projective(system(G), system(F), NullHomogenization(), seed; kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogenous and the other not!")
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end
        homvar = uniquevar(F)
        homogenization = Homogenization(1)
        var_ordering = [homvar; MP.variables(F)]
        Gₕ = system(homogenize(G, homvar), var_ordering)
        Fₕ = system(homogenize(F, homvar), var_ordering)
        Projective(Gₕ, Fₕ, homogenization, seed; kwargs...), prob.startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(prob::ParameterSystem, homvar, seed; system=FPSystem, kwargs...)
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
            error(overdetermined_error_msg)
        end
        if homvar === nothing
            F, variables, NullHomogenization()
        else
            F, variables, Homogenization(homvar, variables)
        end
    else
        # N = n is the only valid size configuration
        if n > N
            error(overdetermined_error_msg)
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
