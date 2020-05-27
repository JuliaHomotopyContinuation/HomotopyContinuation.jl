export AbstractProblem,
    Problem,
    TrackingType,
    AffineTracking,
    ProjectiveTracking,
    homotopy,
    homogenization,
    embed,
    homvars,
    problem_startsolutions


const problem_startsolutions_supported_keywords = [
    [
        :seed,
        :homvar,
        :homvars,
        :variable_groups,
        :homotopy,
        :system,
        :system_scaling,
        :affine_tracking,
        :projective_tracking,
        :start_system,
        :only_torus,
    ]
    input_supported_keywords
]


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

abstract type AbstractProblem{T<:TrackingType,VG<:VariableGroups} end
Base.broadcastable(P::AbstractProblem) = Ref(P)

"""
    homotopy(prob::AbstractProblem)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::AbstractProblem) = prob.homotopy

"""
    seed(prob::AbstractProblem)

Get the random seed used for the problem `prob`.
"""
seed(prob::AbstractProblem) = prob.seed


"""
    Problem(T::TrackingType, H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `Problem`. If `T <: ProjectiveTracking` then the homotopy `H` needs to be homogeneous.
"""
struct Problem{T<:TrackingType,VG<:VariableGroups} <: AbstractProblem{T,VG}
    tracking_type::T
    homotopy::AbstractHomotopy
    vargroups::VG
    seed::Int
    startsolutions_need_reordering::Bool
    regauging_factors::Union{Nothing,Vector{Float64}}
end
function Problem(
    T::TrackingType,
    H::AbstractHomotopy,
    vargroups::VariableGroups,
    seed::Int;
    startsolutions_need_reordering = false,
    regauging_factors = nothing,
)
    if isa(T, AffineTracking) && startsolutions_need_reordering
        error("Affine tracking doesn't support reodering.")
    end
    Problem(T, H, vargroups, seed, startsolutions_need_reordering, regauging_factors)
end
function Problem(
    T::TrackingType,
    G::AbstractSystem,
    F::AbstractSystem,
    vargroups::VariableGroups,
    seed::Int;
    homotopy = DEFAULT_HOMOTOPY,
    kwargs...,
)
    Problem(T, homotopy(G, F), vargroups, seed; kwargs...)
end
Problem{T}(args...; kwargs...) where {T<:TrackingType} = Problem(T(), args...; kwargs...)


struct PolyhedralProblem{T<:TrackingType,VG<:VariableGroups} <: AbstractProblem{T,VG}
    tracking_type::T
    toric_homotopy::AbstractHomotopy
    generic_homotopy::AbstractHomotopy
    vargroups::VG
    seed::Int
    startsolutions_need_reordering::Bool
    regauging_factors::Union{Nothing,Vector{Float64}}
end
function PolyhedralProblem(
    T::TrackingType,
    toric_homotopy::ToricHomotopy,
    generic_homotopy::AbstractHomotopy,
    vargroups::VariableGroups,
    seed::Int;
    startsolutions_need_reordering = false,
    regauging_factors = nothing,
)
    if isa(T, AffineTracking) && startsolutions_need_reordering
        error("Affine tracking doesn't support reodering.")
    end
    PolyhedralProblem(
        T,
        toric_homotopy,
        generic_homotopy,
        vargroups,
        seed,
        startsolutions_need_reordering,
        regauging_factors,
    )
end
PolyhedralProblem{T}(args...; kwargs...) where {T<:TrackingType} =
    Problem(T(), args...; kwargs...)


struct OverdeterminedProblem{
    T<:TrackingType,
    VG<:VariableGroups,
    AP<:AbstractProblem{T,VG},
} <: AbstractProblem{T,VG}
    problem::AP
    target_system::AbstractSystem
    regauging_factors::Union{Nothing,Vector{Float64}}
end
function OverdeterminedProblem(
    problem::AbstractProblem,
    target_system::AbstractSystem;
    regauging_factors = nothing,
)
    OverdeterminedProblem(problem, target_system, regauging_factors)
end


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
homvars(prob::OverdeterminedProblem) = homvars(prob.problem)

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
embed!(x, prob::OverdeterminedProblem{ProjectiveTracking}, v) = embed!(x, prob.problem, v)
embed!(x, prob::AbstractProblem{ProjectiveTracking}, v::PVector) =
    ProjectiveVectors.embed!(x, v)
embed!(x, prob::AbstractProblem{AffineTracking}, v::AbstractVector) = (x .= v; x)

function embed(prob::AbstractProblem{ProjectiveTracking}, v)
    if prob.startsolutions_need_reordering
        embed_projective(prob.vargroups, v)
    else
        ProjectiveVectors.embed(v, projective_dims(prob.vargroups))
    end
end
embed(prob::OverdeterminedProblem{ProjectiveTracking}, v) = embed(prob.problem, v)
embed(prob::AbstractProblem{ProjectiveTracking}, v::PVector) = v
embed(prob::AbstractProblem{AffineTracking}, v::AbstractVector) = v

tracking_vector_type(prob::AbstractProblem{AffineTracking}) = Vector{ComplexF64}
function tracking_vector_type(
    prob::AbstractProblem{ProjectiveTracking,<:VariableGroups{N}},
) where {N}
    PVector{ComplexF64,N}
end

"""
    pull_back(prob::Problem, x)

Pull the solution `x` into affine space if necessary. Creates a copy.
"""
function pull_back(prob::AbstractProblem{AffineTracking}, x::AbstractVector; regauge = true)
    if regauge
        regauge!(copy(x), prob)
    else
        copy(x)
    end
end
function pull_back(prob::OverdeterminedProblem{ProjectiveTracking}, x::PVector; kwargs...)
    pull_back(prob.problem, x; kwargs...)
end
function pull_back(prob::AbstractProblem{ProjectiveTracking}, x::PVector; regauge = true)
    if pull_back_is_to_affine(prob)
        λ = prob.regauging_factors
        if λ === nothing || !regauge
            map(kᵢ -> x[kᵢ[1]] / x[kᵢ[2]], prob.vargroups.pull_back_mapping)
        else
            map(
                kᵢ -> (λ[kᵢ[1]] * x[kᵢ[1]]) / (λ[kᵢ[2]] * x[kᵢ[2]]),
                prob.vargroups.pull_back_mapping,
            )
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
function pull_back_is_to_affine(prob::OverdeterminedProblem{ProjectiveTracking})
    pull_back_is_to_affine(prob.problem)
end

function construct_system(F::Composition, system_constructor; homvars = nothing, kwargs...)
    CompositionSystem(F, system_constructor; homvars = homvars, kwargs...)
end
function construct_system(F::MPPolys, system_constructor; homvars = nothing, kwargs...)
    system_constructor(F; kwargs...)
end

function apply_system_scaling(F, vars, system_scaling::Union{Nothing,Symbol,Bool})
    if system_scaling == :equations_and_variables || system_scaling == true
        precondition(F, vars)
    elseif system_scaling == :equations
        normalize_coefficients(F), nothing
    elseif system_scaling === nothing || system_scaling == false
        F, nothing
    else
        throw(ArgumentError(
            "Got unsupported argument `system_scaling=$(system_scaling)`." *
            " Valid values are `nothing`, `:equations` and `:equations_and_variables`.",
        ))
    end
end

@nospecialize

function default_affine_tracking(F::TargetSystemInput{<:MPPolyInputs}, hominfo)
    !(is_homogeneous(F.system, hominfo))
end
function default_affine_tracking(
    F::TargetSystemInput{<:MPPolyInputs},
    hominfo::HomogenizationInformation,
)
    is_hom = is_homogeneous(F.system, hominfo)
    if !is_hom && !isnothing(hominfo.homvars)
        throw(ArgumentError("Input system is not homogeneous although `homvars=$(hominfo.homvars)` was passed."))
    end

    !(is_hom || (hominfo.vargroups !== nothing && length(hominfo.vargroups) > 1))
end
function default_affine_tracking(input::StartTargetInput{<:MPPolyInputs}, hominfo)
    F, G = input.target, input.start
    !(is_homogeneous(F, hominfo) && is_homogeneous(G, hominfo))
end
function default_affine_tracking(input::HomotopyInput, hominfo)
    !(is_homogeneous(input.H))
end
function default_affine_tracking(input::AbstractHomotopy, hominfo)
    !(is_homogeneous(input))
end
function default_affine_tracking(F::ParameterSystemInput{<:MPPolyInputs}, hominfo)
    !(is_homogeneous(F.system, hominfo; parameters = F.parameters))
end
function default_affine_tracking(input::ParameterSystemInput{<:AbstractSystem}, hominfo)
    H = ParameterHomotopy(
        input.system,
        p₁ = input.p₁,
        p₀ = input.p₀,
        γ₁ = input.γ₁,
        γ₀ = input.γ₀,
    )
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


function problem_startsolutions(
    args...;
    seed = randseed(),
    variable_ordering = nothing,
    kwargs...,
)
    Random.seed!(seed)
    supported, rest = splitkwargs(kwargs, input_supported_keywords)
    input, startsolutions =
        input_startsolutions(args...; variable_ordering = variable_ordering, supported...)
    if variable_ordering !== nothing
        problem_startsolutions(
            input,
            startsolutions,
            seed;
            variable_ordering = variable_ordering,
            rest...,
        )
    else
        problem_startsolutions(input, startsolutions, seed; rest...)
    end
end

function problem_startsolutions(
    input::AbstractInput,
    startsolutions = nothing;
    seed = randseed(),
    kwargs...,
)
    problem_startsolutions(input, startsolutions, seed; kwargs...)
end
function problem_startsolutions(
    input::AbstractInput,
    startsolutions,
    seed::Int;
    homvar = nothing,
    homvars = nothing,
    variable_groups = nothing,
    affine_tracking = nothing,
    projective_tracking = nothing,
    kwargs...,
)


    homvar_info = HomogenizationInformation(;
        homvar = homvar, homvars = homvars, variable_groups = variable_groups,
    )

    #projective_tracking is for the frontend
    #internally, we use affine_tracking
    if isnothing(affine_tracking) && !isnothing(projective_tracking)
        affine_tracking = !projective_tracking
    end

    if isnothing(affine_tracking)
        problem_startsolutions(input, startsolutions, homvar_info, seed; kwargs...)
    else
        problem_startsolutions(
            input,
            startsolutions,
            homvar_info,
            seed;
            affine_tracking = affine_tracking,
            kwargs...,
        )
    end

end

function problem_startsolutions(
    input::HomotopyInput,
    startsolutions,
    homvar_info,
    seed;
    affine_tracking = default_affine_tracking(input, homvar_info),
    system_scaling = nothing,
    variable_ordering = nothing,
    kwargs...,
)
    if !affine_tracking
        is_homogeneous(input.H) || throw(ArgumentError("Input homotopy is not homogeneous by our numerical check."))
        if homvar_info === nothing
            N = size(input.H)[2]
            if start_solution_sample(startsolutions) isa PVector
                vargroups = VariableGroups(N)
            else
                vargroups = VariableGroups(N, N)
            end
            Problem{ProjectiveTracking}(input.H, vargroups, seed), startsolutions
        else
            vargroups = VariableGroups(size(input.H)[2], homvar_info)
            Problem{ProjectiveTracking}(input.H, vargroups, seed), startsolutions
        end
    else
        Problem{AffineTracking}(input.H, VariableGroups(size(input.H)[2]), seed),
        startsolutions
    end
end



##############
# TOTALDEGREE
##############

function problem_startsolutions(
    input::TargetSystemInput{<:MPPolyInputs},
    ::Nothing,
    homvar_info,
    seed;
    start_system = :total_degree,
    affine_tracking = default_affine_tracking(input, homvar_info),
    system_scaling = DEFAULT_SYSTEM_SCALING,
    system = DEFAULT_SYSTEM,
    only_torus = false,
    variable_ordering = variables(input.system),
    kwargs...,
)
    if affine_tracking
        vargroups = VariableGroups(variable_ordering)
        homvars = nothing
        F = input.system
        tracking_type = AffineTracking()
    else
        F, vargroups, homvars =
            homogenize_if_necessary(input.system, homvar_info; vars = variable_ordering)
        tracking_type = ProjectiveTracking()
    end

    classification = classify_system(F, vargroups; affine_tracking = affine_tracking)
    if classification == :underdetermined
        throw(ArgumentError(
            "Underdetermined polynomial systems are currently not supported." *
            " Consider adding linear polynomials to your system in order to reduce your system" *
            " to a zero dimensional system.",
        ))
# The following case is too annoying right now
    elseif classification == :overdetermined && ngroups(vargroups) > 1
        throw(ArgumentError(
            "Overdetermined polynomial systems with a multi-homogeneous" *
            " structure are currently not supported.",
        ))
    end

    vars = flattened_variable_groups(vargroups)
    target_constructor(f) = construct_system(f, system; variables = vars, homvars = homvars)

    # Scale systems
    f, regauging_factors = apply_system_scaling(F, vars, system_scaling)

    if ngroups(vargroups) > 1 # Multihomogeneous
        affine_tracking && error("Affine tracking is currently not supported for variable groups.")

        D = multidegrees(F, vargroups)
        C = multi_bezout_coefficients(D, projective_dims(vargroups))

        g = MultiHomTotalDegreeSystem(D, C)
        problem = Problem(
            tracking_type,
            g,
            target_constructor(f),
            vargroups,
            seed;
            regauging_factors = regauging_factors, kwargs...,
        )
        startsolutions = totaldegree_solutions(g, vargroups)
    # We have ngroups(vargroups) == 1
    elseif start_system == :total_degree
        degrees = maxdegrees(F)

        if classification == :overdetermined
            # We need to square up the system. This will only be a homogenous system if
            # all the degrees are identical
            if !affine_tracking && !all(d -> d == degrees[1], degrees)
                throw(ArgumentError("Cannot square up a homogenous system where not all polynomial systems have a the same degree. Consider starting with an affine system."))
            end
            perm = sortperm(degrees; rev = true)
            # reorder polynomials by minimal degree
            f = permute(f, perm)
            degrees = degrees[perm]
        end
        n = affine_tracking ? length(vars) : length(vars) - 1

        g = TotalDegreeSystem(degrees[1:n]; affine = affine_tracking)

        if classification == :overdetermined
            A = randn(ComplexF64, n, length(f) - n)
            f̄ = target_constructor(f)
            problem = OverdeterminedProblem(
                Problem(
                    tracking_type,
                    g,
                    SquaredUpSystem(f̄, A, degrees),
                    vargroups,
                    seed;
                    kwargs...,
                ),
                f̄;
                regauging_factors = regauging_factors,
            )
        # square case
        else
            problem = Problem(
                tracking_type,
                g,
                target_constructor(f),
                vargroups,
                seed;
                regauging_factors = regauging_factors, kwargs...,
            )
        end
        startsolutions = totaldegree_solutions(g, vargroups)
    elseif start_system == :polyhedral
        n = affine_tracking ? length(vars) : length(vars) - 1
        degrees = maxdegrees(F)
        if f isa Composition
            f = expand(f)
        end

        if classification == :overdetermined
            # square the system up
            # We do not automatically sort the polynomials by degree since it is not
            # clear that this is actually always advantagous
            f̂ = [LinearAlgebra.I randn(ComplexF64, n, length(f) - n)] * f
        else
            f̂ = f
        end

        support, coeffs = support_coefficients(f̂, vars)
        if !affine_tracking
            affine_support = map(A -> A[1:end-1, :], support)
        else
            affine_support = support
        end

        if !only_torus
            for i = 1:length(support)
                !iszero(@view affine_support[i][:, end]) || continue
                affine_support[i] = [affine_support[i] zeros(Int32, n)]
                push!(coeffs[i], zero(eltype(coeffs[i])))
                if !affine_tracking
                    support[i] = [support[i] [zeros(Int32, n); degrees[i]]]
                end
            end
        end
        cell_iter = PolyhedralStartSolutionsIterator(affine_support)
        toric_homotopy = ToricHomotopy(affine_support, cell_iter.start_coefficients)
        generic_homotopy =
            CoefficientHomotopy(support, cell_iter.start_coefficients, coeffs)

        if classification == :overdetermined
            # always use fp system for the target system since
            problem = OverdeterminedProblem(
                PolyhedralProblem(
                    tracking_type,
                    toric_homotopy,
                    generic_homotopy,
                    vargroups,
                    seed;
                    kwargs...,
                ),
                FPSystem(f),
                regauging_factors = regauging_factors,
            )
        else
            problem = PolyhedralProblem(
                tracking_type,
                toric_homotopy,
                generic_homotopy,
                vargroups,
                seed;
                regauging_factors = regauging_factors, kwargs...,
            )
        end
        startsolutions = cell_iter
    else
        throw(ArgumentError(
            "Unsupported argument `start_system=$start_system`. " *
            "Possible values are `:total_degree` and `:polyhedral`",
        ))
    end

    problem, startsolutions
end

@noinline function problem_startsolutions(
    input::TargetSystemInput{ModelKit.System},
    ::Nothing,
    hominfo::Union{Nothing,HomogenizationInformation},
    seed;
    start_system::Symbol = :total_degree,
    variable_ordering = nothing,
    kwargs...,
)
    if start_system === :polyhedral
        a::Tuple{Any,Any} = modelkit_polyhedral(input.system, hominfo, seed; kwargs...)
    else
        problem_startsolutions(
            TargetSystemInput(ModelKitSystem(input.system)),
            nothing,
            hominfo,
            seed;
            start_system = :total_degree,
            variable_ordering = nothing, #variable_ordering,
            kwargs...,
        )
    end
end

@noinline function modelkit_polyhedral(system, hominfo, seed; kwargs...)
    vars = polyvar.(system.variables)
    F = evaluate(system.expressions, system.variables => vars)::Any
    problem_startsolutions(
        TargetSystemInput(F),
        nothing,
        hominfo,
        seed;
        start_system = :polyhedral,
        variable_ordering = vars,
    )::Tuple{Any,Any}
end

@noinline function problem_startsolutions(
    input::TargetSystemInput{<:AbstractSystem},
    ::Nothing,
    hominfo::Union{Nothing,HomogenizationInformation},
    seed;
    start_system = :total_degree,
    system = DEFAULT_SYSTEM,
    affine_tracking = nothing,
    system_scaling = nothing,
    variable_ordering = nothing,
    kwargs...,
)
    start_system == :total_degree || throw(ArgumentError("`start_system = :total_degree` is the only possibility with `AbstractSystem`s."))
    n, N = size(input.system)

    degrees, is_homogeneous = degrees_ishomogeneous(input.system)
    variable_groups = VariableGroups(N, hominfo)
    classification =
        classify_system(input.system, variable_groups; affine_tracking = !is_homogeneous)

    # overdetermined
    if classification == :overdetermined
        if is_homogeneous && !all(d -> d == degrees[1], degrees)
            throw(ArgumentError("Cannot square up a homogenous system where not all polynomial systems have a the same degree. Consider starting with an affine system."))
        end
        if ngroups(variable_groups) > 1
            throw(ArgumentError(
                "Overdetermined polynomial systems with a multi-homogeneous" *
                " structure are currently not supported.",
            ))
        end


        m = N - is_homogeneous
        A = randn(ComplexF64, m, n - m)
        problem = OverdeterminedProblem(
            Problem{is_homogeneous ? ProjectiveTracking : AffineTracking}(
                TotalDegreeSystem(degrees[1:m]; affine = !is_homogeneous),
                SquaredUpSystem(input.system, A, degrees),
                variable_groups,
                seed;
                kwargs...,
            ),
            input.system,
        )
        starts = totaldegree_solutions(max.(degrees[1:m], maximum(degrees[m+1:end])))
    # underdetermined
    elseif classification == :underdetermined
        throw(ArgumentError(
            "Underdetermined polynomial systems are currently not supported." *
            " Consider adding linear polynomials to your system in order to reduce your system" *
            " to a zero dimensional system.",
        ))
    else # square
        G = TotalDegreeSystem(degrees; affine = !is_homogeneous)

        problem = Problem{is_homogeneous ? ProjectiveTracking : AffineTracking}(
            G,
            input.system,
            variable_groups,
            seed;
            kwargs...,
        )
        starts = totaldegree_solutions(degrees)
    end
    # Check overdetermined case
    problem, starts
end

function degrees_ishomogeneous(F)
    n, N = size(F)

    degs = compute_numerically_degrees(F)
    is_homogeneous = !isnothing(degs)
    if isnothing(degs)
        degs = degrees(F)
    end
    isnothing(
        degs,
    ) && throw(ArgumentError("Cannot compute degrees of the input system. Consider overloading `system_degree(F)."))

    degs, is_homogeneous
end

###############
# START TARGET
###############

function problem_startsolutions(
    input::StartTargetInput{<:MPPolyInputs,<:MPPolyInputs},
    startsolutions,
    homvar,
    seed;
    affine_tracking = default_affine_tracking(input, homvar),
    system_scaling = DEFAULT_SYSTEM_SCALING,
    system = DEFAULT_SYSTEM,
    variable_ordering = variables(input.target),
    kwargs...,
)
    vars = variable_ordering
    F, G = input.target, input.start
    F_ishom, G_ishom = is_homogeneous.((F, G))
    if !affine_tracking && F_ishom && G_ishom
        f, regauging_factors = apply_system_scaling(F, vars, system_scaling)
        vargroups = VariableGroups(vars, homvar)
        F̄ = construct_system(f, system; variables = vars, homvars = homvar)
        Ḡ = construct_system(G, system; variables = vars, homvars = homvar)
        prob = Problem{ProjectiveTracking}(
            Ḡ,
            F̄,
            vargroups,
            seed;
            regauging_factors = regauging_factors,
            startsolutions_need_reordering = true,
            kwargs...,
        )
        prob, startsolutions
    elseif F_ishom || G_ishom
        throw(ArgumentError("One of the input polynomials is homogeneous and the other not. Make either both homogeneous or non-homogeneous!"))
    elseif !affine_tracking
        homvar === nothing || throw(ArgumentError("Input system is not homogeneous although `homvar` was passed."))

        h = uniquevar(F)
        push!(vars, h)
        sort!(vars, rev = true)
        g, f = homogenize(G, h), homogenize(F, h)
        f, regauging_factors = apply_system_scaling(f, vars, system_scaling)
        F̄ = construct_system(f, system; variables = vars, homvars = homvar)
        Ḡ = construct_system(g, system; variables = vars, homvars = homvar)
        vargroups = VariableGroups(vars, h)
        prob = Problem{ProjectiveTracking}(
            Ḡ,
            F̄,
            vargroups,
            seed;
            regauging_factors = regauging_factors,
            startsolutions_need_reordering = true,
            kwargs...,
        )
        prob, startsolutions
    else # affine_tracking
        f, regauging_factors = apply_system_scaling(F, vars, system_scaling)
        F̂ = construct_system(f, system; variables = vars)
        Ĝ = construct_system(G, system; variables = vars)
        prob = Problem{AffineTracking}(
            Ĝ,
            F̂,
            VariableGroups(vars),
            seed;
            regauging_factors = regauging_factors, kwargs...,
        )
        prob, startsolutions
    end
end

#####################
# Parameter homotopy
#####################

function problem_startsolutions(
    input::ParameterSystemInput{<:MPPolyInputs},
    startsolutions,
    hominfo,
    seed;
    affine_tracking = default_affine_tracking(input, hominfo),
    system = SPSystem,
    system_scaling = nothing,
    variable_ordering = variables(input.system, input.parameters),
    kwargs...,
)
    Prob = affine_tracking ? Problem{AffineTracking} : Problem{ProjectiveTracking}
    if affine_tracking
        variable_groups = VariableGroups(variable_ordering)
        homvars = nothing
        F = input.system
    else
        F, variable_groups, homvars = homogenize_if_necessary(
            input.system,
            hominfo;
            vars = variable_ordering, parameters = input.parameters,
        )
    end
    vars = flattened_variable_groups(variable_groups)
    F̂ = construct_system(
        F,
        system;
        homvars = homvars, variables = vars, parameters = input.parameters,
    )
    H = ParameterHomotopy(F̂, p₁ = input.p₁, p₀ = input.p₀, γ₁ = input.γ₁, γ₀ = input.γ₀)
    Prob(H, variable_groups, seed; startsolutions_need_reordering = !affine_tracking),
    startsolutions
end


function problem_startsolutions(
    @nospecialize(input::ParameterSystemInput{<:AbstractSystem}),
    startsolutions,
    hominfo,
    seed;
    affine_tracking = default_affine_tracking(input, hominfo),
    system = nothing,
    system_scaling = nothing,
    kwargs...,
)
    n, N = size(input.system)
    H = ParameterHomotopy(
        input.system,
        p₁ = input.p₁,
        p₀ = input.p₀,
        γ₁ = input.γ₁,
        γ₀ = input.γ₀,
    )
    variable_groups = VariableGroups(N, hominfo)
    tracking_type = affine_tracking ? AffineTracking() : ProjectiveTracking()
    Problem(
        tracking_type,
        H,
        variable_groups,
        seed;
        startsolutions_need_reordering = false,
    ),
    startsolutions
end

function problem_startsolutions(
    input::ParameterSystemInput{ModelKit.System},
    startsolutions,
    hominfo,
    seed;
    affine_tracking = nothing,
    system = nothing,
    system_scaling = nothing,
    kwargs...,
)
    n, N = size(input.system)

    # construct parameter homotopy
    has_gamma = input.γ₁ !== nothing && input.γ₀ !== nothing
    params = ComplexF64[input.p₁; input.p₀]
    has_gamma && push!(params, input.γ₁, input.γ₀)
    H = ModelKitHomotopy(
        build_parameter_homotopy(input.system; gamma = has_gamma);
        parameters = params,
    )
    variable_groups = VariableGroups(N, hominfo)
    if affine_tracking === nothing
        affine_tracking = default_affine_tracking(H, hominfo)
    end

    variable_groups = VariableGroups(N, hominfo)
    tracking_type = affine_tracking ? AffineTracking() : ProjectiveTracking()
    Problem(
        tracking_type,
        H,
        variable_groups,
        seed;
        startsolutions_need_reordering = false,
    ),
    startsolutions
end

function build_parameter_homotopy(F::ModelKit.System; gamma = false)
    p = F.parameters
    if gamma
        @var __start_params__[1:length(p)] __target_params__[1:length(p)] __t__ __γ__[1:2]
        a = __γ__[1] * __t__
        b = (1 - __t__) * __γ__[2]
        h = subs(
            F.expressions,
            p => (a .* __start_params__ .+ b .* __target_params__) ./ (a + b),
        )
        Homotopy(h, F.variables, __t__, [__start_params__; __target_params__; __γ__])
    else
        @var __start_params__[1:length(p)] __target_params__[1:length(p)] __t__
        h = subs(
            F.expressions,
            p => __t__ .* __start_params__ .+ (1 - __t__) .* __target_params__,
        )
        Homotopy(h, F.variables, __t__, [__start_params__; __target_params__])
    end
end

@specialize
