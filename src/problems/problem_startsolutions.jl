export problem_startsolutions

const supported_keywords = [[:seed, :homvar, :homotopy, :system, :scale_systems]; Input.supported_keywords]
const DEFAULT_SYSTEM = FPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

function construct_system(F::Composition, system_constructor; homvars=nothing, kwargs...)
	Systems.CompositionSystem(F, system_constructor; homvars=homvars, kwargs...)
end
function construct_system(F::MPPolys, system_constructor; homvars=nothing, kwargs...)
	system_constructor(F; kwargs...)
end

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
function problem_startsolutions(input::AbstractInput, seed;
	homvar::Union{Nothing, Int, MP.AbstractVariable}=nothing, kwargs...)
    problem_startsolutions(input, homvar, seed; kwargs...)
end

function problem_startsolutions(input::Input.Homotopy, homvar, seed; kwargs...)
    Projective(input.H, VariableGroups(size(input.H)[2], homvar), seed), input.startsolutions
end


##############
# TOTALDEGREE
##############

function problem_startsolutions(prob::TotalDegree{<:Input.MPPolyInputs}, homvar, seed; scale_systems=true, system=DEFAULT_SYSTEM, kwargs...)
    F, variables, variable_groups, homvars = Utilities.homogenize_if_necessary(prob.system; homvar=homvar)

	check_square_homogenous_system(F, variable_groups)

	# Scale systems
	if scale_systems
		G = homogenous_totaldegree_polysystem(prob.degrees, variables, variable_groups)
		_, f, G_scaling_factors, _ = Utilities.scale_systems(G, F, report_scaling_factors=true)
		g = Systems.TotalDegreeSystem(prob.degrees, G_scaling_factors)
		problem = Projective(g,
						construct_system(f, system; variables=variables, homvars=homvars),
						variable_groups, seed; kwargs...)
	else
		problem = Projective(Systems.TotalDegreeSystem(prob.degrees),
			construct_system(F, system; variables=variables, homvars=homvars), variable_groups, seed; kwargs...)

	end
	startsolutions = totaldegree_solutions(prob.degrees)

	problem, startsolutions
end

function homogenous_totaldegree_polysystem(degrees, variables, variable_groups::VariableGroups{1})
	vars = variables[variable_groups.groups[1]]
	map(1:length(degrees)) do i
		d = degrees[i]
		vars[i]^d - vars[end]^d
	end
end

function problem_startsolutions(prob::TotalDegree{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
    G = Systems.TotalDegreeSystem(prob.degrees)
	# Check overdetermined case
	n > N && error(Utilities.overdetermined_error_msg)
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

function problem_startsolutions(prob::StartTarget{<:Input.MPPolyInputs, <:Input.MPPolyInputs}, homvar, seed; scale_systems=true, system=DEFAULT_SYSTEM, kwargs...)
    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
	vars = variables(F)
    if F_ishom && G_ishom
		vargroups = VariableGroups(vars, homvar)
		if scale_systems
			g, f = Utilities.scale_systems(G, F)
		else
			g, f = G, F
		end
		F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system; variables=vars, homvars=homvar)
        Projective(Ḡ, F̄, vargroups, seed; kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogenous and the other not!")
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end

		h = uniquevar(F)
        push!(vars, h)
        sort!(vars, rev=true)
		if scale_systems
			g, f = Utilities.scale_systems(homogenize(G, h), homogenize(F, h))
		else
			g, f = homogenize(G, h), homogenize(F, h)
		end
        F̄ = construct_system(f, system, variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system, variables=vars, homvars=homvar)
		vargroups = VariableGroups(vars, h)
        Projective(Ḡ, F̄, vargroups, seed; kwargs...), prob.startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(prob::ParameterSystem, homvar, seed; system=SPSystem, kwargs...)
    F, variables, variable_groups, homvars = Utilities.homogenize_if_necessary(prob.system, homvar=homvar, parameters=prob.parameters)
	F̄ = construct_system(F, system; homvars=homvars, variables=variables, parameters=prob.parameters)
    H = ParameterHomotopy(F̄, p₁=prob.p₁, p₀=prob.p₀, γ₁=prob.γ₁, γ₀=prob.γ₀)

    Projective(H, variable_groups, seed), prob.startsolutions
end
