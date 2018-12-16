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
    F, variables, variable_groups = Utilities.homogenize_if_necessary(prob.system)
	# Since homvar is provided we either need to homogenize or
	# we have already a homogenous system.
    if variable_groups.dedicated_homvars # affine case
		# Check overdetermined case
		length(F) ≥ length(variables) && error(overdetermined_error_msg)

        proj = Projective(
            Systems.TotalDegreeSystem(prob.degrees),
            system(F, variables), variable_groups, seed; kwargs...)
        proj, totaldegree_solutions(prob.degrees, homogenous=false)
    else
		# Check overdetermined case
		length(F) > length(variables) && error(overdetermined_error_msg)

        G = Systems.TotalDegreeSystem(prob.degrees)
        start = totaldegree_solutions(prob.degrees, homogenous=true)
        Projective(G, system(F, variables), variable_groups, seed; kwargs...), start
    end
end

function problem_startsolutions(prob::TotalDegree{Vector{AP}},
    homvar::MP.AbstractVariable, seed; system=DEFAULT_SYSTEM, kwargs...) where {AP<:MP.AbstractPolynomialLike}

    if !ishomogenous(prob.system)
        error("Input system is not homogenous although `homvar=$(homvar)` was passed.")
    end
    F, variables, variable_groups = Utilities.homogenize_if_necessary(prob.system; homvar=homvar)
	# Check overdetermined case
	length(F) > length(variables) && error(overdetermined_error_msg)

    start = totaldegree_solutions(prob.degrees, homogenous=false)
    proj = Projective(
        Systems.TotalDegreeSystem(prob.degrees),
        system(F, variables), variable_groups, seed; kwargs...)
    proj, start
end


function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    G = Systems.TotalDegreeSystem(prob.degrees)
	# Check overdetermined case
	n > N && error(overdetermined_error_msg)
	variable_groups = VariableGroups(N, homvaridx)
    (Projective(G, prob.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(prob.degrees))
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Int, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)

    G = Systems.TotalDegreeSystem(prob.degrees)
	variable_groups = VariableGroups(N, homvaridx)
    (Projective(G, prob.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(prob.degrees))
end


###############
# START TARGET
###############

function problem_startsolutions(prob::StartTarget{Vector{AP1}, Vector{AP2}}, homvar, seed; system=DEFAULT_SYSTEM, kwargs...) where
    {AP1<:MP.AbstractPolynomialLike, AP2<:MP.AbstractPolynomialLike}

    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
	variables = MP.variables(F)
    if F_ishom && G_ishom
		vargroups = VariableGroups(variables, homvar)
        Projective(system(G), system(F), vargroups, seed; kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogenous and the other not!")
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end

		h = uniquevar(F)
        push!(variables, h)
        sort!(variables, rev=true)
        F′ = homogenize(F, h)
		G′ = homogenize(G, h)

		vargroups = VariableGroups(variables, h)
        Projective(system(G′, variables), system(F′, variables), vargroups, seed; kwargs...), prob.startsolutions
    end
end
#
# #####################
# # Parameter homotopy
# #####################
#
# function problem_startsolutions(prob::ParameterSystem, homvar, seed; system=FPSystem, kwargs...)
#     F, variables, homogenization = Utilities.homogenize_if_necessary(prob.system, homvar=homvar, parameters=prob.parameters)
#
#     H = ParameterHomotopy(F, prob.parameters, variables=variables,
# 						  p₁=prob.p₁, p₀=prob.p₀, γ₁=prob.γ₁, γ₀=prob.γ₀)
#
#     Projective(H, homogenization, seed), prob.startsolutions
# end
