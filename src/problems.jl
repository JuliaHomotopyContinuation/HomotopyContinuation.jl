export AbstractProblem, ProjectiveProblem, homotopy, homogenization, embed, homvars,
	problem_startsolutions


const problem_startsolutions_supported_keywords = [
	[:seed, :homvar, :homvars, :variable_groups, :homotopy, :system, :system_scaling, :affine];
	input_supported_keywords]


const DEFAULT_SYSTEM = @static VERSION < v"1.1" ? FPSystem : SPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

abstract type AbstractProblem end
Base.broadcastable(P::AbstractProblem) = Ref(P)
"""
    homotopy(prob::AbstractProblem)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::AbstractProblem) = prob.homotopy

"""
    ProjectiveProblem(H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `ProjectiveProblemProblem`. The homotopy `H` needs to be homogenous.
"""
struct ProjectiveProblem{H<:AbstractHomotopy, VG<:VariableGroups} <: AbstractProblem
    homotopy::H
    vargroups::VG
    seed::Int
	startsolutions_need_reordering::Bool
end
function ProjectiveProblem(H::AbstractHomotopy, vargroups::VariableGroups, seed::Int; startsolutions_need_reordering=false)
	ProjectiveProblem(H, vargroups, seed, startsolutions_need_reordering)
end
function ProjectiveProblem(G::AbstractSystem, F::AbstractSystem,
        vargroups::VariableGroups, seed::Int;
		homotopy=DEFAULT_HOMOTOPY, kwargs...)
    ProjectiveProblem(homotopy(G, F), vargroups, seed; kwargs...)
end

"""
    homvars(prob::ProjectiveProblem)

Get the homogenization variables of the problem. Returns `nothing` if there are no.
"""
function homvars(prob::ProjectiveProblem)
    if has_dedicated_homvars(prob.vargroups)
        map(last, prob.vargroups.groups)
    else
        nothing
    end
end


"""
    AffineProblem(H::AbstractHomotopy, seed::Int)

Construct a `AffineProblem`.
"""
struct AffineProblem{H<:AbstractHomotopy, T} <: AbstractProblem
    homotopy::H
    vargroups::VariableGroups{1, false, T}
    seed::Int
end
function AffineProblem(G::AbstractSystem, F::AbstractSystem,
        vargroups::VariableGroups{1, false}, seed::Int; homotopy=DEFAULT_HOMOTOPY)
    AffineProblem(homotopy(G, F), vargroups, seed)
end


"""
    embed(prob::ProjectiveProblem, v)

Embed the vector `v` into projective space if necessary.
"""
function embed!(x, prob::ProjectiveProblem, v)
	if prob.startsolutions_need_reordering
		embed_projective!(x, prob.vargroups, v)
	else
		ProjectiveVectors.embed!(x, v)
	end
end
embed!(x, prob::ProjectiveProblem, v::PVector) = ProjectiveVectors.embed!(x, v)
embed!(x, prob::AffineProblem, v::AbstractVector) = begin x .= v; x end

function embed(prob::ProjectiveProblem, v)
	if prob.startsolutions_need_reordering
		embed_projective(prob.vargroups, v)
	else
		ProjectiveVectors.embed(v, projective_dims(prob.vargroups))
	end
end
embed(prob::ProjectiveProblem, v::PVector) = v
embed(prob::AffineProblem, v::AbstractVector) = v

"""
    pull_back(prob::ProjectiveProblem, x)

Pull the solution `x` into affine space if necessary. Creates a copy.
"""
pull_back(prob::ProjectiveProblem, x::PVector) = pull_back(prob.vargroups, x)
pull_back(prob::AffineProblem, x::AbstractVector) = copy(x)

function construct_system(F::Composition, system_constructor; homvars=nothing, kwargs...)
	CompositionSystem(F, system_constructor; homvars=homvars, kwargs...)
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
    * `system=SPSystem`: A constructor to assemble a [`AbstractSystem`](@ref). The constructor
    is called with `system(polys, variables)` where `variables` determines the variable ordering.
    * `homotopy=StraightLineHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref) an `AbstractSystem`. The constructor
    is called with `homotopy(start, target)` where `start` and `target` are systems constructed
    with `system`.
"""
function problem_startsolutions end

function problem_startsolutions(args...; seed=randseed(), kwargs...)
    Random.seed!(seed)
    supported, rest = splitkwargs(kwargs, input_supported_keywords)
    problem_startsolutions(input(args...; supported...), seed; rest...)
end

function problem_startsolutions(input::AbstractInput; seed=randseed(), kwargs...)
    problem_startsolutions(input, seed; kwargs...)
end
function problem_startsolutions(input::AbstractInput, seed;
	homvar=nothing, homvars=nothing, variable_groups=nothing, kwargs...)
	homvar_info = HomogenizationInformation(;homvar=homvar, homvars=homvars, variable_groups=variable_groups)
    problem_startsolutions(input, homvar_info, seed; kwargs...)
end

function problem_startsolutions(input::HomotopyInput, homvar, seed; kwargs...)
    ProjectiveProblem(input.H, VariableGroups(size(input.H)[2], homvar), seed), input.startsolutions
end



##############
# TOTALDEGREE
##############

function problem_startsolutions(prob::TotalDegreeInput{<:MPPolyInputs}, homvar_info, seed; system_scaling=true, system=DEFAULT_SYSTEM, kwargs...)
    F, vargroups, homvars = homogenize_if_necessary(prob.system, homvar_info)
	variables = flattened_variable_groups(vargroups)
	target_constructor = f -> begin
		construct_system(f, system; variables=variables, homvars=homvars)
	end

	check_square_homogenous_system(F, vargroups)

	if ngroups(vargroups) == 1
		degrees = maxdegrees(F)
		# Scale systems
		if system_scaling
			G = homogenous_totaldegree_polysystem(degrees, variables, vargroups)
			_, f, G_scaling_factors, _ = scale_systems(G, F, report_scaling_factors=true)
			g = TotalDegreeSystem(degrees, G_scaling_factors)
			problem = ProjectiveProblem(g, target_constructor(f), vargroups, seed; kwargs...)
		else
			g = TotalDegreeSystem(degrees)
			problem = ProjectiveProblem(g, target_constructor(F), vargroups, seed; kwargs...)
		end
	# Multihomogenous
	else
		D = multidegrees(F, vargroups)
		C = multi_bezout_coefficients(D, projective_dims(vargroups))
		if system_scaling
			G = totaldegree_polysystem(D, vargroups, C)
			_, f, G_scaling_factors, _ = scale_systems(G, F, report_scaling_factors=true)
			g = MultiHomTotalDegreeSystem(D, C, G_scaling_factors)
			problem = ProjectiveProblem(g, target_constructor(f), vargroups, seed; kwargs...)
		else
			g = MultiHomTotalDegreeSystem(D, C)
			problem = ProjectiveProblem(g, target_constructor(F), vargroups, seed; kwargs...)
		end

		problem = ProjectiveProblem(g, target_constructor(F), vargroups, seed; kwargs...)
	end
	startsolutions = totaldegree_solutions(g, vargroups)

	problem, startsolutions
end

function homogenous_totaldegree_polysystem(degrees, variables, variable_groups::VariableGroups{1})
	vars = variables[variable_groups.groups[1]]
	map(1:length(degrees)) do i
		d = degrees[i]
		vars[i]^d - vars[end]^d
	end
end

function problem_startsolutions(prob::TotalDegreeInput{<:AbstractSystem}, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
	degrees = abstract_system_degrees(prob.system)
    G = TotalDegreeSystem(degrees)
	# Check overdetermined case
	n > N && error(overdetermined_error_msg)
	variable_groups = VariableGroups(N, homvaridx)
    (ProjectiveProblem(G, prob.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function problem_startsolutions(prob::TotalDegreeInput{<:AbstractSystem}, hominfo::HomogenizationInformation, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(prob.system)
	degrees = abstract_system_degrees(prob.system)
    G = TotalDegreeSystem(degrees)
	variable_groups = VariableGroups(N, hominfo)
    (ProjectiveProblem(G, prob.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function abstract_system_degrees(F)
	n, N = size(F)
	degrees = check_homogenous_degrees(F)
	# system needs to be homogenous
	if n + 1 > N
		error(overdetermined_error_msg)
	elseif  n + 1 ≠ N
		error("Input system is not a square homogenous system!")
	end
	degrees
end

###############
# START TARGET
###############

function problem_startsolutions(prob::StartTargetInput, homvar, seed; system_scaling=true, system=DEFAULT_SYSTEM, kwargs...)
    F, G = prob.target, prob.start
    F_ishom, G_ishom = ishomogenous.((F, G))
	vars = variables(F)
    if F_ishom && G_ishom
		vargroups = VariableGroups(vars, homvar)
		if system_scaling
			g, f = scale_systems(G, F)
		else
			g, f = G, F
		end
		F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system; variables=vars, homvars=homvar)
        ProjectiveProblem(Ḡ, F̄, vargroups, seed; startsolutions_need_reordering=true, kwargs...), prob.startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogenous and the other not!")
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end

		h = uniquevar(F)
        push!(vars, h)
        sort!(vars, rev=true)
		if system_scaling
			g, f = scale_systems(homogenize(G, h), homogenize(F, h))
		else
			g, f = homogenize(G, h), homogenize(F, h)
		end
        F̄ = construct_system(f, system, variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system, variables=vars, homvars=homvar)
		vargroups = VariableGroups(vars, h)
        ProjectiveProblem(Ḡ, F̄, vargroups, seed; startsolutions_need_reordering=true, kwargs...), prob.startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(prob::ParameterSystemInput, hominfo, seed; affine=false, system=SPSystem, kwargs...)
	if affine
		vars = variables(prob.system; parameters=prob.parameters)
		variable_groups = VariableGroups(vars)
		F = construct_system(prob.system, system; variables=vars, parameters=prob.parameters)
	    H = ParameterHomotopy(F, p₁=prob.p₁, p₀=prob.p₀, γ₁=prob.γ₁, γ₀=prob.γ₀)

	    AffineProblem(H, variable_groups, seed), prob.startsolutions
	else
	    F, variable_groups, homvars = homogenize_if_necessary(prob.system, hominfo; parameters=prob.parameters)
		vars = flattened_variable_groups(variable_groups)
		F̄ = construct_system(F, system; homvars=homvars, variables=vars, parameters=prob.parameters)
	    H = ParameterHomotopy(F̄, p₁=prob.p₁, p₀=prob.p₀, γ₁=prob.γ₁, γ₀=prob.γ₀)

	    ProjectiveProblem(H, variable_groups, seed; startsolutions_need_reordering=true), prob.startsolutions
	end
end
