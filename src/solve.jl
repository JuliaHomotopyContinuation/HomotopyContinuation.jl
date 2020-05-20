export solve, Solver, solver_startsolutions, paths_to_track

struct Solver{T<:AbstractPathTracker}
    tracker::T
    seed::Union{Nothing,UInt32}
end
Solver(tracker::AbstractPathTracker) = Solver(tracker, nothing)

function solver_startsolutions(
    F::Union{System,AbstractSystem},
    starts = nothing;
    seed = rand(UInt32),
    start_system = isnothing(variable_groups(F)) ? :polyhedral : :total_degree,
    p₁ = nothing,
    start_parameters = p₁,
    p₀ = nothing,
    target_parameters = p₀,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)

    if start_parameters !== nothing
        tracker = parameter_homotopy(
            F;
            start_parameters = start_parameters,
            target_parameters = target_parameters,
            kwargs...,
        )
    elseif start_system == :polyhedral
        tracker, starts = polyhedral(F; target_parameters = target_parameters, kwargs...)
    elseif start_system == :total_degree
        tracker, starts = total_degree(F; target_parameters = target_parameters, kwargs...)
    else
        throw(KeywordArgumentException(
            :start_system,
            start_system,
            "Possible values are: `:polyhedral` and `:total_degree`.",
        ))
    end

    Solver(tracker, seed), starts
end

function solver_startsolutions(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem},
    starts = nothing;
    seed = rand(UInt32),
    tracker_options = TrackerOptions(),
    path_tracker_options = PathTrackerOptions(),
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)
    H = start_target_homotopy(G, F; kwargs...)
    tracker =
        PathTracker(H; tracker_options = tracker_options, options = path_tracker_options)

    Solver(tracker, seed), starts
end

function parameter_homotopy(
    F::Union{System,AbstractSystem};
    p₁ = nothing,
    start_parameters = p₁,
    p₀ = nothing,
    target_parameters = p₀,
    tracker_options = TrackerOptions(),
    path_tracker_options = PathTrackerOptions(),
)
    isnothing(start_parameters) && throw(UndefKeywordError(:start_parameters))
    isnothing(target_parameters) && throw(UndefKeywordError(:target_parameters))
    m, n = size(F)
    H = ParameterHomotopy(F, start_parameters, target_parameters)
    f = System(F)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups,) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    PathTracker(H; tracker_options = tracker_options, options = path_tracker_options)
end

function start_target_homotopy(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem};
    γ = 1.0,
    gamma = γ,
)
    f, g = System(F), System(G)


    size(F) == size(G) || error("The provided systems don't have the same size.")
    is_homogeneous(f) == is_homogeneous(g) ||
        error("The provided systems are not both homogeneous.")
    variable_groups(f) == variable_groups(g) ||
        error("The provided systems don't decalare the same variable groups.")

    m, n = size(F)

    H = StraightLineHomotopy(G, F; gamma = gamma)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups,) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    H
end

"""
    solve(...)

TODO
"""
function solve(args...; kwargs...)
    solver, starts = solver_startsolutions(args...; kwargs...)
    solve(solver, starts)
end


function solve(S::Solver, starts)
    path_results = PathResult[]
    for (k, s) in enumerate(starts)
        push!(path_results, track(S.tracker, s; path_number = k))
    end

    Result(path_results; seed = S.seed)
end
(solver::Solver)(starts) = solve(solver, starts)


"""
    paths_to_track(
        f::Union{System,AbstractSystem};
        start_system::Symbol = :polyhedral,
        kwargs...)

Returns the number of paths tracked when calling [`solve`](@ref) with the given arguments.
"""
function paths_to_track(
    f::Union{System,AbstractSystem};
    start_system::Symbol = :polyhedral,
    kwargs...,
)
    paths_to_track(f, Val(start_system); kwargs...)
end
