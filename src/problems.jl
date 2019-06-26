export AbstractProblem, Problem, TrackingType, AffineTracking, ProjectiveTracking,
	   homotopy, homogenization, embed, homvars, problem_startsolutions


const problem_startsolutions_supported_keywords = [
	[:seed, :homvar, :homvars, :variable_groups, :homotopy, :system, :system_scaling,
	:affine_tracking, :projective_tracking, :start_system, :only_torus];
	input_supported_keywords]


const DEFAULT_SYSTEM = @static VERSION < v"1.1" ? FPSystem : SPSystem
const DEFAULT_HOMOTOPY = StraightLineHomotopy
const DEFAULT_SYSTEM_SCALING = :equations

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

abstract type AbstractProblem{T<:TrackingType} end
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
struct Problem{T<:TrackingType, VG<:VariableGroups} <: AbstractProblem{T}
	tracking_type::T
    homotopy::AbstractHomotopy
    vargroups::VG
    seed::Int
	startsolutions_need_reordering::Bool
	regauging_factors::Union{Nothing,Vector{Float64}}
end
function Problem(T::TrackingType, H::AbstractHomotopy, vargroups::VariableGroups, seed::Int;
		startsolutions_need_reordering=false, regauging_factors=nothing)
	if isa(T, AffineTracking) && startsolutions_need_reordering
		error("Affine tracking doesn't support reodering.")
	end
	Problem(T, H, vargroups, seed, startsolutions_need_reordering, regauging_factors)
end
function Problem(T::TrackingType, G::AbstractSystem, F::AbstractSystem,
        vargroups::VariableGroups, seed::Int;
		homotopy=DEFAULT_HOMOTOPY, kwargs...)
    Problem(T, homotopy(G, F), vargroups, seed; kwargs...)
end
Problem{T}(args...; kwargs...) where {T<:TrackingType} = Problem(T(), args...; kwargs...)

"""
    Problem(T::TrackingType, H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `Problem`. If `T <: ProjectiveTracking` then the homotopy `H` needs to be homogeneous.
"""
struct PolyhedralProblem{T<:TrackingType, VG<:VariableGroups} <: AbstractProblem{T}
	tracking_type::T
    toric_homotopy::AbstractHomotopy
	generic_homotopy::AbstractHomotopy
    vargroups::VG
    seed::Int
	startsolutions_need_reordering::Bool
	regauging_factors::Union{Nothing,Vector{Float64}}
end
function PolyhedralProblem(T::TrackingType,
			toric_homotopy::ToricHomotopy,
			generic_homotopy::AbstractHomotopy,
			vargroups::VariableGroups, seed::Int;
		startsolutions_need_reordering=false, regauging_factors=nothing)
	if isa(T, AffineTracking) && startsolutions_need_reordering
		error("Affine tracking doesn't support reodering.")
	end
	PolyhedralProblem(T, toric_homotopy, generic_homotopy,
				vargroups, seed, startsolutions_need_reordering, regauging_factors)
end
PolyhedralProblem{T}(args...; kwargs...) where {T<:TrackingType} = Problem(T(), args...; kwargs...)

"""
    homvars(prob::AbstractProblem)

Get the homogenization variables of the problem. Returns `nothing` if there are no.
"""
function homvars(prob::AbstractProblem{ProjectiveTracking})
    if has_dedicated_homvars(prob.vargroups)
        map(last, prob.vargroups.groups)
    else
        nothing
    end
end
homvars(prob::AbstractProblem{AffineTracking}) = nothing

"""
    embed(prob::Problem{ProjectiveTracking}, v)

Embed the vector `v` into projective space if necessary.
"""
function embed!(x, prob::AbstractProblem{ProjectiveTracking}, v)
	if prob.startsolutions_need_reordering
		embed_projective!(x, prob.vargroups, v)
	else
		ProjectiveVectors.embed!(x, v)
	end
end
embed!(x, prob::AbstractProblem{ProjectiveTracking}, v::PVector) = ProjectiveVectors.embed!(x, v)
embed!(x, prob::AbstractProblem{AffineTracking}, v::AbstractVector) = begin x .= v; x end

function embed(prob::AbstractProblem{ProjectiveTracking}, v)
	if prob.startsolutions_need_reordering
		embed_projective(prob.vargroups, v)
	else
		ProjectiveVectors.embed(v, projective_dims(prob.vargroups))
	end
end
embed(prob::AbstractProblem{ProjectiveTracking}, v::PVector) = v
embed(prob::AbstractProblem{AffineTracking}, v::AbstractVector) = v

"""
    pull_back(prob::Problem, x)

Pull the solution `x` into affine space if necessary. Creates a copy.
"""
function pull_back(prob::AbstractProblem{AffineTracking}, x::AbstractVector; regauge=true)
	if regauge
		regauge!(copy(x), prob)
	else
		copy(x)
	end
end
function pull_back(prob::AbstractProblem{ProjectiveTracking}, x::PVector; regauge=true)
	if pull_back_is_to_affine(prob)
		λ = prob.regauging_factors
		if λ === nothing || !regauge
			map(kᵢ -> x[kᵢ[1]] / x[kᵢ[2]],
				prob.vargroups.pull_back_mapping)
		else
			map(kᵢ -> (λ[kᵢ[1]] * x[kᵢ[1]]) / (λ[kᵢ[2]] * x[kᵢ[2]]),
				prob.vargroups.pull_back_mapping)
		end
	elseif regauge
		LinearAlgebra.normalize!(regauge!(copy(x), prob))
	else
		LinearAlgebra.normalize(x)
	end
end

function regauge!(x, prob::AbstractProblem)
	if prob.regauging_factors !== nothing
		x .= prob.regauging_factors .* x
	end
	x
end

"""
    pull_back_is_to_affine(prob::ProjectiveProblem)

Returns `true` if [`pull_back`](@ref) would pull a solution `x` into affine space.
"""
pull_back_is_to_affine(prob::AbstractProblem{AffineTracking}) = true
function pull_back_is_to_affine(prob::AbstractProblem{ProjectiveTracking})
	has_dedicated_homvars(prob.vargroups)
end

function construct_system(F::Composition, system_constructor; homvars=nothing, kwargs...)
	CompositionSystem(F, system_constructor; homvars=homvars, kwargs...)
end
function construct_system(F::MPPolys, system_constructor; homvars=nothing, kwargs...)
	system_constructor(F; kwargs...)
end

function apply_system_scaling(F, vars, system_scaling::Union{Nothing, Symbol, Bool})
	if system_scaling == :equations_and_variables || system_scaling == true
		precondition(F, vars)
	elseif system_scaling == :equations
		normalize_coefficients(F), nothing
	elseif system_scaling === nothing || system_scaling == false
		F, nothing
	else
		throw(ArgumentError("Got unsupported argument `system_scaling=$(system_scaling)`." *
		  " Valid values are `nothing`, `:equations` and `:equations_and_variables`."))
	end
end

function default_affine_tracking(F::TargetSystemInput{<:MPPolyInputs}, hominfo)
	!(is_homogeneous(F.system, hominfo))
end
function default_affine_tracking(F::TargetSystemInput{<:MPPolyInputs}, hominfo::HomogenizationInformation)
	is_hom = is_homogeneous(F.system, hominfo)
	if !is_hom && !isnothing(hominfo.homvars)
        throw(ArgumentError("Input system is not homogeneous although `homvars=$(hominfo.homvars)` was passed."))
	end

	!(is_hom ||
	  (hominfo.vargroups !== nothing && length(hominfo.vargroups) > 1))
end
function default_affine_tracking(input::StartTargetInput{<:MPPolyInputs}, hominfo)
	F, G = input.target, input.start
	!(is_homogeneous(F, hominfo) && is_homogeneous(G, hominfo))
end
function default_affine_tracking(input::HomotopyInput, hominfo)
	!(is_homogeneous(input.H))
end
function default_affine_tracking(F::ParameterSystemInput{<:MPPolyInputs}, hominfo)
	!(is_homogeneous(F.system, hominfo; parameters=F.parameters))
end
function default_affine_tracking(input::ParameterSystemInput{<:AbstractSystem}, hominfo)
	H = ParameterHomotopy(input.system, p₁=input.p₁, p₀=input.p₀, γ₁=input.γ₁, γ₀=input.γ₀)
	!(is_homogeneous(H))
end


is_homogeneous(F::AbstractSystem) = compute_numerically_degrees(F) !== nothing
is_homogeneous(H::AbstractHomotopy) = is_homogeneous(FixedHomotopy(H, rand()))

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
	* `affine_tracking::Bool=true`: Indicate whether path tracking should happen in affine space.
	* `projective_tracking::Bool=false`: Indicate whether path tracking should happen in projective space. The flag `affine_tracking` is dominant.
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
	homvar=nothing, homvars=nothing, variable_groups=nothing, affine_tracking=nothing, projective_tracking=nothing, kwargs...)

	homvar_info = HomogenizationInformation(;homvar=homvar, homvars=homvars, variable_groups=variable_groups)

	#projective_tracking is for the frontend
	#internally, we use affine_tracking
	if isnothing(affine_tracking)
		if isnothing(projective_tracking)
			problem_startsolutions(input, startsolutions, homvar_info, seed; kwargs...)
		else
			problem_startsolutions(input, startsolutions, homvar_info, seed; affine_tracking=!projective_tracking, kwargs...)
		end
	else
		problem_startsolutions(input, startsolutions, homvar_info, seed; affine_tracking=affine_tracking, kwargs...)
	end
end

function problem_startsolutions(input::HomotopyInput, startsolutions, homvar_info, seed; affine_tracking=default_affine_tracking(input, homvar_info), system_scaling=nothing, kwargs...)
	if !affine_tracking
		is_homogeneous(input.H) || throw(ArgumentError("Input homotopy is not homogeneous by our numerical check."))
		if homvar_info === nothing
			N = size(input.H)[2]
			if start_solution_sample(startsolutions) isa PVector
				vargroups = VariableGroups(N)
			else
				vargroups = VariableGroups(N,N)
			end
    		Problem{ProjectiveTracking}(input.H, vargroups, seed), startsolutions
		else
			vargroups = VariableGroups(size(input.H)[2], homvar_info)
			Problem{ProjectiveTracking}(input.H, vargroups, seed), startsolutions
		end
	else
		Problem{AffineTracking}(input.H, VariableGroups(size(input.H)[2]), seed), startsolutions
	end
end



##############
# TOTALDEGREE
##############

function problem_startsolutions(input::TargetSystemInput{<:MPPolyInputs}, ::Nothing, homvar_info, seed;
				start_system=:total_degree,
				affine_tracking=default_affine_tracking(input, homvar_info),
				system_scaling=DEFAULT_SYSTEM_SCALING, system=DEFAULT_SYSTEM,
				only_torus=false, kwargs...)
	if affine_tracking
		vargroups = VariableGroups(variables(input.system))
		homvars = nothing
		F = input.system
		tracking_type = AffineTracking()
	else
		F, vargroups, homvars = homogenize_if_necessary(input.system, homvar_info)
		tracking_type = ProjectiveTracking()
	end

	classifcation = classify_system(F, vargroups; affine_tracking=affine_tracking)
	if classifcation == :underdetermined
		throw(ArgumentError("Underdetermined polynomial systems are currently not supported." *
		     " Consider adding linear polynomials to your system in order to reduce your system" *
			 " to a zero dimensional system."))
	# The following case is too annoying right now
	elseif classifcation == :overdetermined && ngroups(vargroups) > 1
		error("Overdetermined polynomial systems with a multi-homogenous structure are currently not supported.")
    end

	vars = flattened_variable_groups(vargroups)
	target_constructor(f) = construct_system(f, system; variables=vars, homvars=homvars)

	# Scale systems
	f, regauging_factors = apply_system_scaling(F, vars, system_scaling)

	if ngroups(vargroups) > 1 # Multihomogeneous
		affine_tracking && error("Affine tracking is currently not supported for variable groups.")

		D = multidegrees(F, vargroups)
		C = multi_bezout_coefficients(D, projective_dims(vargroups))

		g = MultiHomTotalDegreeSystem(D, C)
		problem = Problem(tracking_type, g, target_constructor(f), vargroups, seed;
					regauging_factors=regauging_factors, kwargs...)
		startsolutions = totaldegree_solutions(g, vargroups)
    # We have ngroups(vargroups) == 1
	elseif start_system == :total_degree
		n = affine_tracking ? length(vars) : length(vars) - 1
		degrees = maxdegrees(F)

		if classifcation == :overdetermined
			perm = sortperm(degrees; rev=true)
			# reorder polynomials by minimal degree
			F = permute(F, perm)
			degrees = degrees[perm]
		end

		if classifcation == :overdetermined
			A = randn(ComplexF64, n, npolynomials(F) - n)
			f̂ = SquaredUpSystem(target_constructor(f), A, degrees)
		else # square case
			f̂ = target_constructor(f)
		end
		g = TotalDegreeSystem(degrees[1:n]; affine=affine_tracking)
		problem = Problem(tracking_type, g, f̂, vargroups, seed;
					regauging_factors=regauging_factors, kwargs...)
		startsolutions = totaldegree_solutions(g, vargroups)
	elseif start_system == :polyhedral
		n = affine_tracking ? length(vars) : length(vars) - 1
		degrees = maxdegrees(F)
		if classifcation == :overdetermined
			error("Overdetermined systems are only supported with `algorithm=:total_degree`")
		end
		support, coeffs = support_coefficients(f, vars)
		if !affine_tracking
			affine_support = map(A -> A[1:end-1,:], support)
		else
			affine_support = support
		end

		if !only_torus
			for i in 1:length(support)
				if !iszero(@view affine_support[i][:,end])
					affine_support[i] = [affine_support[i] zeros(Int32, n)]
					push!(coeffs[i], zero(eltype(coeffs[i])))
					if !affine_tracking
						support[i] = [support[i] [zeros(Int32, n); degrees[i]]]
					end
				end
			end
		end

		cell_iter = PolyhedralStartSolutionsIterator(affine_support)
		toric_homotopy = ToricHomotopy(affine_support, cell_iter.start_coefficients)
	    generic_homotopy = CoefficientHomotopy(support, cell_iter.start_coefficients, coeffs)

		problem = PolyhedralProblem(tracking_type, toric_homotopy, generic_homotopy, vargroups, seed;
					regauging_factors=regauging_factors, kwargs...)
		startsolutions = cell_iter
	else
		throw(ArgumentError("Unsupported argument `start_system=$start_system`. " *
            "Possible values are `:total_degree` and `:polyhedral`"))
	end

	problem, startsolutions
end

function problem_startsolutions(input::TargetSystemInput{<:AbstractSystem},
				::Nothing, homvaridx::Nothing, seed;
				system=DEFAULT_SYSTEM, system_scaling=nothing, kwargs...)
    n, N = size(input.system)

	degrees = abstract_system_degrees(input.system)
    G = TotalDegreeSystem(degrees)
	# Check overdetermined case
	n > N && throw(ArgumentError("Currently custom overdetermined systems are not supported."))
	variable_groups = VariableGroups(N, homvaridx)
    (Problem{ProjectiveTracking}(G, input.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function problem_startsolutions(input::TargetSystemInput{<:AbstractSystem},
				::Nothing, hominfo::HomogenizationInformation, seed;
				system=DEFAULT_SYSTEM, system_scaling=nothing, kwargs...)
    n, N = size(input.system)
	degrees = abstract_system_degrees(input.system)
    G = TotalDegreeSystem(degrees)
	variable_groups = VariableGroups(N, hominfo)
    (Problem{ProjectiveTracking}(G, input.system, variable_groups, seed; kwargs...),
     totaldegree_solutions(degrees))
end

function abstract_system_degrees(F)
	n, N = size(F)

	degrees = compute_numerically_degrees(F)
	isnothing(degrees) && throw(ArgumentError("Input system is not homogeneous by our numerical check."))
	# system needs to be homogeneous
	if n + 1 > N
		throw(ArgumentError("Currently custom overdetermined systems are not supported."))
	elseif  n + 1 ≠ N
		throw(ArgumentError("Input system is not a square homogeneous system!"))
	end
	degrees
end

###############
# START TARGET
###############

function problem_startsolutions(input::StartTargetInput{<:MPPolyInputs, <:MPPolyInputs}, startsolutions, homvar, seed;
				affine_tracking=default_affine_tracking(input, homvar),
				system_scaling=DEFAULT_SYSTEM_SCALING,
				system=DEFAULT_SYSTEM, kwargs...)
    F, G = input.target, input.start
    F_ishom, G_ishom = is_homogeneous.((F, G))
	vars = variables(F)
    if !affine_tracking && F_ishom && G_ishom
		f, regauging_factors = apply_system_scaling(F, vars, system_scaling)
		vargroups = VariableGroups(vars, homvar)
		F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(G, system; variables=vars, homvars=homvar)
        prob = Problem{ProjectiveTracking}(Ḡ, F̄, vargroups, seed;
				regauging_factors=regauging_factors,
				startsolutions_need_reordering=true, kwargs...)
		prob, startsolutions
    elseif F_ishom || G_ishom
        throw(ArgumentError("One of the input polynomials is homogeneous and the other not. Make either both homogeneous or non-homogeneous!"))
    elseif !affine_tracking
        homvar === nothing || throw(ArgumentError("Input system is not homogeneous although `homvar` was passed."))

		h = uniquevar(F)
        push!(vars, h)
        sort!(vars, rev=true)
		g, f = homogenize(G, h), homogenize(F, h)
		f, regauging_factors = apply_system_scaling(f, vars, system_scaling)
        F̄ = construct_system(f, system; variables=vars, homvars=homvar)
		Ḡ = construct_system(g, system; variables=vars, homvars=homvar)
		vargroups = VariableGroups(vars, h)
        prob = Problem{ProjectiveTracking}(Ḡ, F̄, vargroups, seed;
				regauging_factors=regauging_factors,
				startsolutions_need_reordering=true, kwargs...)
		prob, startsolutions
	else # affine_tracking
		f, regauging_factors = apply_system_scaling(F, vars, system_scaling)
		F̂ = construct_system(f, system; variables=vars)
		Ĝ = construct_system(G, system; variables=vars)
		prob = Problem{AffineTracking}(Ĝ, F̂, VariableGroups(vars), seed;
				regauging_factors=regauging_factors, kwargs...)
		prob, startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(input::ParameterSystemInput{<:MPPolyInputs}, startsolutions, hominfo, seed;
						affine_tracking=default_affine_tracking(input, hominfo),
						system=SPSystem, system_scaling=nothing, kwargs...)
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

function problem_startsolutions(input::ParameterSystemInput{<:AbstractSystem}, startsolutions, hominfo, seed;
				affine_tracking=default_affine_tracking(input, hominfo),
				system=nothing, system_scaling=nothing, kwargs...)
	n, N = size(input.system)
	H = ParameterHomotopy(input.system, p₁=input.p₁, p₀=input.p₀, γ₁=input.γ₁, γ₀=input.γ₀)
	variable_groups = VariableGroups(N, hominfo)
	tracking_type = affine_tracking ? AffineTracking() : ProjectiveTracking()
	Problem(tracking_type, H, variable_groups, seed; startsolutions_need_reordering=false), startsolutions
end
