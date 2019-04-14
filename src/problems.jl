export AbstractProblem, Problem, TrackingType, AffineTracking, ProjectiveTracking,
	   homotopy, homogenization, embed, homvars, problem_startsolutions


const problem_startsolutions_supported_keywords = [
	[:seed, :homvar, :homvars, :variable_groups, :homotopy, :system, :system_scaling, :affine_tracking];
	input_supported_keywords]


const DEFAULT_SYSTEM = @static VERSION < v"1.1" ? FPSystem : SPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy

"""
	TrackingType

Abstract type determining how paths are tracked.
"""
abstract type TrackingType end

"""
	AffineTracking <: TrackingType

Indicates that paths should be tracked in affine space.
"""
struct AffineTracking <: TrackingType end

"""
	ProjectiveTracking <: TrackingType

Indicates that paths should be tracked in projective space.
"""
struct ProjectiveTracking <: TrackingType end

abstract type AbstractProblem end
Base.broadcastable(P::AbstractProblem) = Ref(P)
"""
    homotopy(prob::AbstractProblem)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::AbstractProblem) = prob.homotopy

"""
    Problem(T::TrackingType, H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `Problem`. If `T <: ProjectiveTracking` then the homotopy `H` needs to be homogeneous.
"""
struct Problem{T<:TrackingType, VG<:VariableGroups} <: AbstractProblem
	tracking_type::T
    homotopy::AbstractHomotopy
    vargroups::VG
    seed::Int
	startsolutions_need_reordering::Bool
end
function Problem(T::TrackingType, H::AbstractHomotopy, vargroups::VariableGroups, seed::Int; startsolutions_need_reordering=false)
	if isa(T, Problem{AffineTracking}) && startsolutions_need_reordering
		error("Affine tracking doesn't support reodering.")
	end
	Problem(T, H, vargroups, seed, startsolutions_need_reordering)
end
function Problem(T::TrackingType, G::AbstractSystem, F::AbstractSystem,
        vargroups::VariableGroups, seed::Int;
		homotopy=DEFAULT_HOMOTOPY, kwargs...)
    Problem(T, homotopy(G, F), vargroups, seed; kwargs...)
end
Problem{T}(args...; kwargs...) where {T<:TrackingType} = Problem(T(), args...; kwargs...)


"""
    homvars(prob::AbstractProblem)

Get the homogenization variables of the problem. Returns `nothing` if there are no.
"""
function homvars(prob::Problem{ProjectiveTracking})
    if has_dedicated_homvars(prob.vargroups)
        map(last, prob.vargroups.groups)
    else
        nothing
    end
end
homvars(prob::Problem{AffineTracking}) = nothing

"""
    embed(prob::Problem{ProjectiveTracking}, v)

Embed the vector `v` into projective space if necessary.
"""
function embed!(x, prob::Problem{ProjectiveTracking}, v)
	if prob.startsolutions_need_reordering
		embed_projective!(x, prob.vargroups, v)
	else
		ProjectiveVectors.embed!(x, v)
	end
end
embed!(x, prob::Problem{ProjectiveTracking}, v::PVector) = ProjectiveVectors.embed!(x, v)
embed!(x, prob::Problem{AffineTracking}, v::AbstractVector) = begin x .= v; x end

function embed(prob::Problem{ProjectiveTracking}, v)
	if prob.startsolutions_need_reordering
		embed_projective(prob.vargroups, v)
	else
		ProjectiveVectors.embed(v, projective_dims(prob.vargroups))
	end
end
embed(prob::Problem{ProjectiveTracking}, v::PVector) = v
embed(prob::Problem{AffineTracking}, v::AbstractVector) = v

"""
    pull_back(prob::Problem, x)

Pull the solution `x` into affine space if necessary. Creates a copy.
"""
pull_back(prob::Problem{AffineTracking}, x::AbstractVector) = copy(x)
pull_back(prob::Problem{ProjectiveTracking}, x::PVector) = pull_back(prob.vargroups, x)
function pull_back(VG::VariableGroups{M, false}, v::PVector{<:Number, M}) where {M}
	LinearAlgebra.normalize(v)
end
function pull_back(VG::VariableGroups{M, true}, x::PVector{<:Number, M}) where {M}
	map(ki -> x[ki[1]] / x[ki[2]], VG.pull_back_mapping)
end

"""
    pull_back_is_to_affine(prob::ProjectiveProblem, x)

Returns `true` if [`pull_back`](@ref) would pull the solution `x` into affine space.
"""
pull_back_is_to_affine(prob::Problem{AffineTracking}, x::AbstractVector) = true
function pull_back_is_to_affine(prob::Problem{ProjectiveTracking}, x::PVector)
	pull_back_is_to_affine(prob.vargroups, x)
end
pull_back_is_to_affine(::VariableGroups{M,true}, ::PVector{<:Number, M}) where {M} = true
pull_back_is_to_affine(::VariableGroups{M,false},::PVector{<:Number, M}) where {M} = false


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
	input, startsolutions = input_startsolutions(args...; supported...)
	problem_startsolutions(input, startsolutions, seed; rest...)
end

function problem_startsolutions(input::AbstractInput, startsolutions=nothing; seed=randseed(), kwargs...)
    problem_startsolutions(input, startsolutions, seed; kwargs...)
end
function problem_startsolutions(input::AbstractInput, startsolutions, seed::Int;
	homvar=nothing, homvars=nothing, variable_groups=nothing, kwargs...)
	homvar_info = HomogenizationInformation(;homvar=homvar, homvars=homvars, variable_groups=variable_groups)
    problem_startsolutions(input, startsolutions, homvar_info, seed; kwargs...)
end

function problem_startsolutions(input::HomotopyInput, startsolutions, homvar, seed; affine_tracking=false, kwargs...)
	if !affine_tracking
		check_homogeneous_degrees(FixedHomotopy(input.H, rand()))
    	Problem{ProjectiveTracking}(input.H, VariableGroups(size(input.H)[2], homvar), seed), startsolutions
	else
		Problem{AffineTracking}(input.H, VariableGroups(size(input.H)[2]), seed), startsolutions
	end
end



##############
# TOTALDEGREE
##############

function problem_startsolutions(input::TotalDegreeInput{<:MPPolyInputs}, ::Nothing, homvar_info, seed; affine_tracking=false, system_scaling=true, system=DEFAULT_SYSTEM, kwargs...)
	if affine_tracking
		vargroups = VariableGroups(variables(input.system))
		homvars = nothing
		F = input.system
		Prob = Problem{AffineTracking}
	else
		F, vargroups, homvars = homogenize_if_necessary(input.system, homvar_info)
		Prob = Problem{ProjectiveTracking}
	end

	classifcation = classify_system(F, vargroups; affine_tracking=affine_tracking)
	if classifcation == :underdetermined
        error("Underdetermined polynomial systems are currently not supported.")
	# The following case is too annoying right now
	elseif classifcation == :overdetermined && ngroups(vargroups) > 1
		error("Overdetermined polynomial systems with a multi-homogenous structure are currently not supported.")
    end

	vars = flattened_variable_groups(vargroups)
	target_constructor(f) = construct_system(f, system; variables=vars, homvars=homvars)

	if ngroups(vargroups) == 1
		n = affine_tracking ? length(vars) : length(vars) - 1
		degrees = maxdegrees(F)

		if classifcation == :overdetermined
			perm = sortperm(degrees; rev=true)
			# reorder polynomials by minimal degree
			F = F[perm]
			degrees = degrees[perm]
		end

		# Scale systems
		if system_scaling && classifcation == :square
			G = totaldegree_polysystem(degrees, vars, vargroups; affine_tracking=affine_tracking)
			_, f, G_scaling_factors, _ = scale_systems(G, F, report_scaling_factors=true)
			g = TotalDegreeSystem(degrees, G_scaling_factors; affine=affine_tracking)
		elseif classifcation == :overdetermined
			g, f = TotalDegreeSystem(degrees[1:n]; affine=affine_tracking), F
		else
			g, f = TotalDegreeSystem(degrees; affine=affine_tracking), F
		end
		Prob = affine_tracking ? Problem{AffineTracking} : Problem{ProjectiveTracking}

		if classifcation == :overdetermined
			A = randn(ComplexF64, n, length(F) - n)
			f̂ = SquaredUpSystem(target_constructor(f), A, degrees)
		else # square case
			f̂ = target_constructor(f)
		end

		problem = Prob(g, f̂, vargroups, seed; kwargs...)
	# Multihomogeneous
	else
		affine_tracking && error("Affine tracking is currently not supported for variable groups.")

		D = multidegrees(F, vargroups)
		C = multi_bezout_coefficients(D, projective_dims(vargroups))
		if system_scaling
			G = totaldegree_polysystem(D, vargroups, C)
			_, f, G_scaling_factors, _ = scale_systems(G, F, report_scaling_factors=true)
			g = MultiHomTotalDegreeSystem(D, C, G_scaling_factors)
			problem = Problem{ProjectiveTracking}(g, target_constructor(f), vargroups, seed; kwargs...)
		else
			g = MultiHomTotalDegreeSystem(D, C)
			problem = Problem{ProjectiveTracking}(g, target_constructor(F), vargroups, seed; kwargs...)
		end

		problem = Problem{ProjectiveTracking}(g, target_constructor(F), vargroups, seed; kwargs...)
	end
	startsolutions = totaldegree_solutions(g, vargroups)

	problem, startsolutions
end

function totaldegree_polysystem(degrees, variables, variable_groups::VariableGroups{1}; affine_tracking=false)
	vars = variables[variable_groups.groups[1]]
	map(1:length(degrees)) do i
		d = degrees[i]
		affine_tracking ? vars[i]^d - 1 : vars[i]^d - vars[end]^d
	end
end

function problem_startsolutions(input::TotalDegreeInput{<:AbstractSystem}, ::Nothing, homvaridx::Nothing, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(input.system)
	degrees = abstract_system_degrees(input.system)
    G = TotalDegreeSystem(degrees)
	# Check overdetermined case
	n > N && error(overdetermined_error_msg)
	variable_groups = VariableGroups(N, homvaridx)
    (Problem{ProjectiveTracking}(G, input.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function problem_startsolutions(input::TotalDegreeInput{<:AbstractSystem}, ::Nothing, hominfo::HomogenizationInformation, seed; system=DEFAULT_SYSTEM, kwargs...)
    n, N = size(input.system)
	degrees = abstract_system_degrees(input.system)
    G = TotalDegreeSystem(degrees)
	variable_groups = VariableGroups(N, hominfo)
    (Problem{ProjectiveTracking}(G, input.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function abstract_system_degrees(F)
	n, N = size(F)
	degrees = check_homogeneous_degrees(F)
	# system needs to be homogeneous
	if n + 1 > N
		error(overdetermined_error_msg)
	elseif  n + 1 ≠ N
		error("Input system is not a square homogeneous system!")
	end
	degrees
end

###############
# START TARGET
###############

function problem_startsolutions(input::StartTargetInput, startsolutions, homvar, seed; affine_tracking=false, system_scaling=true, system=DEFAULT_SYSTEM, kwargs...)
    F, G = input.target, input.start
    F_ishom, G_ishom = ishomogeneous.((F, G))
	vars = variables(F)
    if !affine_tracking && F_ishom && G_ishom
		vargroups = VariableGroups(vars, homvar)
		if system_scaling
			g, f = scale_systems(G, F)
		else
			g, f = G, F
		end
		F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system; variables=vars, homvars=homvar)
        Problem{ProjectiveTracking}(Ḡ, F̄, vargroups, seed; startsolutions_need_reordering=true, kwargs...), startsolutions
    elseif F_ishom || G_ishom
        error("One of the input polynomials is homogeneous and the other not!")
	elseif affine_tracking && F_ishom && G_ishom
		error("Both of the input systems are homogenous but `affine_tracking=true` was passed.")
    elseif !affine_tracking
        homvar === nothing || error("Input system is not homogeneous although `homvar` was passed.")

		h = uniquevar(F)
        push!(vars, h)
        sort!(vars, rev=true)
		if system_scaling
			g, f = scale_systems(homogenize(G, h), homogenize(F, h))
		else
			g, f = homogenize(G, h), homogenize(F, h)
		end
        F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system; variables=vars, homvars=homvar)
		vargroups = VariableGroups(vars, h)
        Problem{ProjectiveTracking}(Ḡ, F̄, vargroups, seed; startsolutions_need_reordering=true, kwargs...), startsolutions
	else # affine_tracking
		F̂ = construct_system(F, system; variables=vars)
		Ĝ = construct_system(G, system; variables=vars)
		Problem{AffineTracking}(Ĝ, F̂, VariableGroups(vars), seed; kwargs...), startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(input::ParameterSystemInput{<:MPPolyInputs}, startsolutions, hominfo, seed; affine_tracking=false, system=SPSystem, kwargs...)
	Prob = affine_tracking ? Problem{AffineTracking} : Problem{ProjectiveTracking}
	if affine_tracking
		variable_groups = VariableGroups(variables(input.system; parameters=input.parameters))
		homvars = nothing
		F = input.system
	else
		F, variable_groups, homvars = homogenize_if_necessary(input.system, hominfo; parameters=input.parameters)
	end
	vars = flattened_variable_groups(variable_groups)
	F̂ = construct_system(F, system; homvars=homvars, variables=vars, parameters=input.parameters)
	H = ParameterHomotopy(F̂, p₁=input.p₁, p₀=input.p₀, γ₁=input.γ₁, γ₀=input.γ₀)
	Prob(H, variable_groups, seed; startsolutions_need_reordering=!affine_tracking), startsolutions
end

function problem_startsolutions(input::ParameterSystemInput{<:AbstractSystem}, startsolutions, hominfo, seed; affine_tracking=nothing, system=SPSystem, kwargs...)
	n, N = size(input.system)
	H = ParameterHomotopy(input.system, p₁=input.p₁, p₀=input.p₀, γ₁=input.γ₁, γ₀=input.γ₀)
	variable_groups = VariableGroups(N, hominfo)

	# We do not set affine tracking here, to handle the case that the input system
	# is not homogenous.
	if affine_tracking === nothing
	   affine_tracking = compute_numerically_degrees(FixedHomotopy(H, rand())) === nothing
   end

	if affine_tracking
		Problem{AffineTracking}(H, variable_groups, seed; startsolutions_need_reordering=false), startsolutions
	else
    	Problem{ProjectiveTracking}(H, variable_groups, seed; startsolutions_need_reordering=false), startsolutions
	end
end
