export NumericalIrreducibleDecomposition,
    nid,
    numerical_irreducible_decomposition,
    regeneration,
    witness_sets,
    decompose,
    ncomponents,
    n_components,
    seed

"""
    WitnessPoints

A wrapper for storing data about witness points.
"""
mutable struct WitnessPoints{Sub1<:AbstractSubspace,Sub2<:AbstractSubspace}
    L::Sub1
    Lᵤ::Sub2
    R::Vector{Vector{ComplexF64}}
end
linear_subspace(W::WP) where {WP<:WitnessPoints} = W.L
linear_subspace_u(W::WP) where {WP<:WitnessPoints} = W.Lᵤ
points(W::WP) where {WP<:WitnessPoints} = W.R
codim(W::WP) where {WP<:WitnessPoints} = dim(linear_subspace(W))
dim(W::WP) where {WP<:WitnessPoints} = codim(linear_subspace(W))
ModelKit.degree(W::WP) where {WP<:WitnessPoints} = length(points(W))
push!(W::WP, R::Vector{ComplexF64}) where {WP<:WitnessPoints} = push!(W.R, R)

function u_transform(W::WitnessPoints)
    L = linear_subspace_u(W)
    A = extrinsic(L).A
    b = extrinsic(L).b

    E = ExtrinsicDescription(A[:, 1:(end-1)], b; orthonormal = false)
    P = map(points(W)) do p
        p[1:(end-1)]
    end

    P, LinearSubspace(E)
end


Base.@kwdef mutable struct WitnessSetsProgress
    progress_meter::PM.ProgressUnknown
    ambient_dim::Int
    degrees::Dict{Int,Int} = Dict{Int,Int}()
    is_solving::Bool = false
    is_computing_hypersurfaces::Bool = true
    is_membership_test::Bool = false
    is_monodromy::Bool = false
    is_finished::Bool = false
    current_task::Int = 0
    ntasks::Int = 0
    current_hypersurface::Int = 1
    nhypersurfaces::Int
end
WitnessSetsProgress(n::Int, nhypersurfaces::Int, progress_meter::PM.ProgressUnknown) =
    WitnessSetsProgress(
        ambient_dim = n,
        nhypersurfaces = nhypersurfaces,
        progress_meter = progress_meter,
    )
update_progress!(
    progress::Nothing;
    is_computing_hypersurfaces::Union{Nothing,Bool} = nothing,
    is_solving::Union{Nothing,Bool} = nothing,
    is_membership_test::Union{Nothing,Bool} = nothing,
    is_monodromy::Union{Nothing,Bool} = nothing,
) = nothing
function update_progress!(
    progress::WitnessSetsProgress;
    is_computing_hypersurfaces::Union{Nothing,Bool} = nothing,
    is_solving::Union{Nothing,Bool} = nothing,
    is_membership_test::Union{Nothing,Bool} = nothing,
    is_monodromy::Union{Nothing,Bool} = nothing,
)
    isnothing(is_computing_hypersurfaces) ? nothing :
    progress.is_computing_hypersurfaces = is_computing_hypersurfaces
    isnothing(is_solving) ? nothing : progress.is_solving = is_solving
    isnothing(is_membership_test) ? nothing :
    progress.is_membership_test = is_membership_test
    isnothing(is_monodromy) ? nothing : progress.is_monodromy = is_monodromy
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress);
        force = true,
    )
end
update_progress_hypersurface!(progress::Nothing, i::Int) = nothing
function update_progress_hypersurface!(progress::WitnessSetsProgress, i::Int)
    progress.current_hypersurface = i
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
# update_progress_tasks!(progress::Nothing, i::Int, m::Int) = nothing: not needed because already defined
function update_progress_tasks!(progress::WitnessSetsProgress, i::Int, m::Int)
    progress.current_task = i
    progress.ntasks = m
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress!(progress::Nothing, W::Union{WitnessPoints,WitnessSet}) = nothing
function update_progress!(progress::WitnessSetsProgress, W::Union{WitnessPoints,WitnessSet})
    progress.degrees[codim(W)] = degree(W)
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress!(progress::Union{Nothing,WitnessSetsProgress}, W::Nothing) = nothing
finish_progress!(progress::Nothing) = nothing
function finish_progress!(progress::WitnessSetsProgress)
    progress.is_computing_hypersurfaces = false
    progress.is_solving = false
    progress.is_membership_test = false
    progress.is_monodromy = false
    progress.is_finished = true
    PM.finish!(progress.progress_meter, showvalues = showvalues(progress))
end
function showvalues(progress::WitnessSetsProgress)

    if !progress.is_finished
        text = [("Status", "")]
        if progress.is_computing_hypersurfaces
            push!(text, ("Computing hypersurfaces", " "))
        else
            push!(
                text,
                (
                    "Intersect with hypersurface",
                    "$(progress.current_hypersurface) / $(progress.nhypersurfaces)",
                ),
            )
            if progress.is_solving
                push!(text, ("Track paths", "$(progress.current_task)/$(progress.ntasks)"))
            elseif progress.is_monodromy
                push!(
                    text,
                    (
                        "Track paths",
                        "$(progress.ntasks)/$(progress.ntasks) (fill up points...)",
                    ),
                )
            elseif progress.is_membership_test
                push!(
                    text,
                    (
                        "Membership test",
                        "$(progress.current_task)/$(progress.ntasks) points checked",
                    ),
                )
            end
        end
    else
        text = [("Status", "done")]
    end

    push!(text, ("Current number of witness points", " "))
    for c = 1:progress.current_hypersurface
        d = progress.ambient_dim - c
        deg = get(progress.degrees, c, nothing)
        if !isnothing(deg) && deg > 0
            push!(text, ("Dimension $d", "$(deg)"))
        end
    end


    text
end

"""
    RegenerationCache

A cache for [`regeneration`](@ref).
"""
mutable struct RegenerationCache{Sys<:AbstractSystem}
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    x0::Vector{ComplexF64}
    y0::Vector{ComplexF64}
    y::Vector{ComplexF64}

    Fᵢ::Sys
    u::Variable
    i::Int
    codim::Int

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{WitnessSetsProgress,Nothing}
end

function RegenerationCache(u, F, codim, EO, TO, progress)
    m, N = size(F)
    A = zeros(ComplexF64, N-1, N)
    b = zeros(ComplexF64, N-1)
    x0 = zeros(ComplexF64, N)
    y0 = zeros(ComplexF64, m)
    y = zeros(ComplexF64, m)

    RegenerationCache(A, b, x0, y0, y, F, u, 0, codim, EO, TO, progress)
end
function update_Fᵢ!(cache, Fᵢ)
    cache.Fᵢ = Fᵢ
    cache.y0 = zeros(ComplexF64, length(Fᵢ))
    cache.y = zeros(ComplexF64, length(Fᵢ))

    nothing
end
update_i!(cache, i) = cache.i = i

"""
    regeneration(F::System; options...) 

This solves ``F=0`` equation-by-equations and returns a [`WitnessSet`](@ref) for every dimension without decomposing them into irreducible components (witness sets that are not decomposed are also called witness supersets).

The implementation is based on the algorithm [u-regeneration](https://arxiv.org/abs/2206.02869) by Duff, Leykin and Rodriguez. 

### Options

* `sorted = true`: if `true`, the polynomials in F will be sorted by degree in decreasing order. 
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `show_progress = true`: indicate whether a progress bar should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) should be displayed.
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).`
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)
* `seed`: choose the random seed.

### Example

The following example computes witness sets for a union of two circles.

```julia-repl
julia> @var x y
julia> f = (x^2 + y^2 - 1) * ((x-1)^2 + (y-1)^2 - 1)
julia> W = regeneration([f])
1-element Vector{WitnessSet}:
    Witness set for dimension 1 of degree 4  
```         
"""
regeneration(F::System; kwargs...) = _regeneration(deepcopy(F); kwargs...)
regeneration(F::Vector{Expression}; kwargs...) = regeneration(System(F); kwargs...)
function _regeneration(
    F::System;
    sorted::Bool = true,
    max_codim::Union{Int,Nothing} = nothing,
    show_progress::Bool = true,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_cond = 1e12,
    ),
    monodromy_options::MonodromyOptions = MonodromyOptions(;
        trace_test = true,
        parameter_sampler = weighted_normal,
    ),
    show_monodromy_progress::Bool = false,
    threading::Bool = Threads.nthreads() > 1,
    seed = nothing,
    atol = 1e-14,
    rtol = sqrt(eps()),
)

    !isnothing(seed) && Random.seed!(seed)

    # the algorithm is u-regeneration as proposed 
    # by Duff, Leykin and Rodriguez in https://arxiv.org/abs/2206.02869

    vars = variables(F)
    n = size(F, 2) # ambient dimension
    c = size(F, 1) # we can have witness sets of codimesion at most min(c,n)
    expected_max_codim = min(c, n)
    if !isnothing(max_codim) && max_codim < expected_max_codim
        # if max_codim is smaller than the expected codimension we must compute witness points for one more codimension, so that we can remove spurious points
        codim = max_codim + 1
    else
        codim = expected_max_codim
    end
    # u-regeneration adds another variable u to F
    @unique_var u
    push!(vars, u)

    # progress bar
    if show_progress
        progress = WitnessSetsProgress(
            n,
            c,
            PM.ProgressUnknown(
                dt = 0.001,
                desc = "Computing witness sets...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end

    # prepare codim witness sets for the output
    # internally we represent a witness superset by WitnessPoints
    # the i-th witness superset out[i] is for codimension i
    out = initialize_witness_sets(codim, n)

    # we compute witness (super)sets for the hypersurfaces f[1]=0,...,f[c]=0.
    # it is covenient to use the WitnessSet wrapper here, because this also keeps track of the equation
    # as a linear subspace we take the linear subspace for out[1], that sets u=0.
    update_progress!(progress; is_computing_hypersurfaces = true)
    H = initialize_hypersurfaces(F, vars, linear_subspace(out[1]); threading = threading)
    if any(H .== nothing)
        return nothing
    end

    # sort expressions by degree
    if sorted
        σ = sortperm(H, by = ModelKit.degree, rev = true)
        f = expressions(F)[σ]
        H = H[σ]
    else
        f = expressions(F)
    end

    # Initialize a cache
    Fᵢ = fixed(System(f[1:1], variables = vars), compile = false)
    cache = RegenerationCache(u, Fᵢ, codim, endgame_options, tracker_options, progress)

    # now comes the core loop of the algorithm.
    # we start with the first hypersurface f[1]=0 and take its witness superset H[1]
    # then, we for i in 2:c we iteratively intersect all current witness sets with f[i]=0
    update_progress!(progress; is_computing_hypersurfaces = false, is_solving = true)
    begin
        for i = 1:c
            update_progress_hypersurface!(progress, i)
            if i == 1
                # the first step: take the witness superset H[1] as initial witness superset
                begin
                    P = solutions(H[1])
                    W = out[1]
                    for p in P
                        push!(W, p)
                    end
                    update_progress!(progress, W)
                end
            else
                begin
                    update_i!(cache, i)

                    # for all W in out we intersect W with H[i]
                    intersect_all!(
                        out,
                        H,
                        cache;
                        threading = threading,
                        atol = atol,
                        rtol = rtol,
                    )

                    # update Fᵢ
                    Fᵢ = fixed(System(f[1:i], variables = vars), compile = false)
                    update_Fᵢ!(cache, Fᵢ)

                    # running one monodromy per witness set to fill up point
                    fill_up!(
                        out,
                        monodromy_options,
                        cache,
                        show_monodromy_progress,
                        threading,
                    )

                    update_progress!(progress)
                end
            end
        end
    end
    pop!(vars)

    if !isnothing(max_codim) && max_codim < expected_max_codim
        # If max_codim < expected_max_codim we have computed one additional set of witness points. Here we remove them again.
        pop!(out)
    end

    filter!(W -> degree(W) > 0, out)
    update_progress!(progress; is_solving = false, is_monodromy = false)
    if !isempty(out)
        witness_sets = map(out) do W
            update_progress!(progress)
            P, L = u_transform(W)
            WitnessSet(fixed(F, compile = false), L, P)
        end
    else
        return witness_sets = Vector{WitnessSet}()
    end


    finish_progress!(progress)

    return witness_sets
end


function get_flag(iter, L₀) 
    # this gives the flag of linear spaces containing L₀ = {Ax=b} and {Ax=b, u=c} indexed by iter; for i in iter, this returns the linear space obtained by deleting the first (i-1) rows from A and b.
    A₀ = extrinsic(L₀).A # per constructor, A₀ has orthonormal rows
    b₀ = extrinsic(L₀).b
    n = size(A₀, 1) + 1
    m = size(A₀, 2)

    # u regenerates operates on 2 types of linear equations Ax=b
    # type 1 does not use u, so we set the last column to zero
    Aᵤ = [A₀ zeros(n - 1)]
    bᵤ = b₀
    # type 2 sets the linear equation u=c. We add this as the first equation
    c = randn(ComplexF64)
    A = [zeros(1, m) 1.0; A₀ zeros(n - 1)]
    b = [c; b₀]

    map(iter) do i 
        j = i + 1
        # type 1 linear equation for out[i] takes the last i rows of Aᵤ
        Eᵤ = ExtrinsicDescription(Aᵤ[i:end, :], bᵤ[i:end]; orthonormal = true)
        Lᵤ = LinearSubspace(Eᵤ)
        # type 2 linear equation for out[i] takes the last i-1 rows and the first row of A
        E = ExtrinsicDescription(A[[1; j:n], :], b[[1; j:n]]; orthonormal = true)
        L = LinearSubspace(E)
        L, Lᵤ
    end
end
function initialize_witness_sets(codim, n)
    # we need a flag containing a linear space of dimension 1
    L₀ = rand_subspace(n; dim = 1)
    flag = get_flag(1:codim, L₀)

    map(flag) do F 
        L, Lᵤ = F
        WitnessPoints(L, Lᵤ, Vector{Vector{ComplexF64}}())
    end
end
function initialize_hypersurfaces(
    F::System,
    vars,
    L;
    threading::Bool = Threads.nthreads() > 1,
)
    f = expressions(F)
    c = length(f)
    out = Vector{WitnessSet}(undef, c)
    for i = 1:c
        fᵢ = f[i]
        h = fixed(System([fᵢ], variables = vars), compile = false)
        pᵢ, qᵢ = get_num_den(fᵢ) # fᵢ = pᵢ / qᵢ
        res = solve(
            System([pᵢ], variables = vars),
            target_subspace = L;
            start_system = :total_degree,
            show_progress = false,
            threading = threading,
        )
        # if fᵢ is rational we check which zeros of pᵢ are also zeros of fᵢ
        if qᵢ != 1
            res = solve(
                h,
                solutions(res);
                start_subspace = L,
                target_subspace = L,
                show_progress = false,
                threading = threading,
            )
        end
        if isnothing(res)
            return nothing
        end
        S = solutions(res, only_nonsingular = true)
        out[i] = WitnessSet(h, L, S)
    end
    out
end

function intersect_all!(out, H, cache; kwargs...)

    i = cache.i
    codim = cache.codim
    progress = cache.progress

    # the i-th hypersurface
    Hᵢ = H[i]

    # we enumerate reversely, so that we can add points to witness sets that we have already
    # taken care of; i.e., if we intersect Wₖ∩Hᵢ below in the intersect_with_hypersurface function,
    # we have already intersected Wₖ₊₁ ∩ Hᵢ.

    E = enumerate(out)
    for (k, Wₖ) in reverse(E) # k = codim(W) for W in Ws
        if k < i
            if k < codim
                Wₖ₊₁ = out[k+1]
            else
                Wₖ₊₁ = nothing
            end
            # here is the intersection step
            # we add points that do not belong to Wₖ∩Hᵢ to Wₖ₊₁.
            intersect_with_hypersurface!(Wₖ, Hᵢ, Wₖ₊₁, cache; kwargs...)
            update_progress!(progress, Wₖ₊₁)
            update_progress!(progress)
        end
        update_progress!(progress, Wₖ)
    end
end

function fill_up!(out, monodromy_options, cache, show_monodromy_progress, threading)
    progress = cache.progress
    Fᵢ = cache.Fᵢ

    update_progress!(progress; is_membership_test = false, is_solving = false, is_monodromy = true)
    for W in out
        if !isnothing(W) && dim(W) > 0 && degree(W) > 0
            res = monodromy_solve(
                Fᵢ,
                W.R,
                linear_subspace(W);
                options = regeneration_monodromy_options(monodromy_options, W),
                show_progress = show_monodromy_progress,
                threading = threading,
                warning = false,
            )
            # check if solutions remain, and if yes, take unique_points
            if nsolutions(res) == 0
                W.R = Vector{Vector{ComplexF64}}()
            else
                W.R = unique_points(solutions(res))
            end
            update_progress!(progress, W)
        end
    end
end

function regeneration_monodromy_options(M::MonodromyOptions, W)

    MonodromyOptions(;
        permutations = false,
        trace_test = true,
        single_loop_per_start_solution = M.single_loop_per_start_solution,
        check_startsolutions = M.check_startsolutions,
        group_actions = M.group_actions,
        loop_finished_callback = M.loop_finished_callback,
        parameter_sampler = M.parameter_sampler,
        equivalence_classes = M.equivalence_classes,
        trace_test_tol = M.trace_test_tol,
        target_solutions_count = Int(floor(1.5 * degree(W))), # in case a singular solution slips through
        timeout = M.timeout,
        min_solutions = M.min_solutions,
        max_loops_no_progress = M.max_loops_no_progress,
        reuse_loops = M.reuse_loops,
        distance = M.distance,
        unique_points_atol = M.unique_points_atol,
        unique_points_rtol = M.unique_points_rtol,
    )
end

"""
    intersect_with_hypersurface!

This is the core routine of the regeneration algorithm. It intersects a set of [`WitnessPoints`](@ref) with a hypersurface.
"""
function intersect_with_hypersurface!(
    W,
    H,
    X,
    cache;
    threading::Bool = Threads.nthreads() > 1,
    kwargs...,
)
    F = cache.Fᵢ
    progress = cache.progress
    u = cache.u

    P = points(W)
    f = (System(F).expressions) # equations for W
    G = system(H) # H is the hypersurface
    g = first(expressions(System(G))) # equations for H

    # Step 1:
    # we check which points of W are also contained in H
    # the rest is removed from P = points(W) and added to P_next
    # for further processing
    m = .!(is_contained(W, H, cache; threading = threading, kwargs...))
    P_next = manage_initial_points!(P, m, W, progress)
    if isempty(P_next)
        return nothing
    end
    # if isnothing(X), then we should not add points to X
    if isnothing(X)
        return nothing
    end

    # Step 2:
    # the points in P_next are used as starting points for a homotopy.
    # where u^d-1 (u is the extra variable in u-regeneration) is deformed into g 
    Hom, d = set_up_u_homotopy(H, u, W, X, f, g, variables(F))

    tracker = EndgameTracker(
        Hom;
        tracker_options = cache.tracker_options,
        options = cache.endgame_options,
    )

    # the start solutions are the Cartesian product between P_next and the d-th roots of unity.
    start = Iterators.product(P_next, [exp(2 * pi * im * k / d) for k = 0:(d-1)])

    # here comes the loop for tracking
    update_progress_tasks!(progress, 0, length(start))
    update_progress!(progress; is_membership_test = false, is_solving = true, is_monodromy = false)
    if threading
        threaded_intersection!(X, start, tracker, progress)
    else
        serial_intersection!(X, start, tracker, progress)
    end
    nothing
end

function manage_initial_points!(P, m, W, progress)
    P_next = P[m]
    deleteat!(P, m)
    return P_next
end

function serial_intersection!(X, start, egtracker, progress)
    l_start = length(start)

    for (i, s) in enumerate(start)
        p = s[1]
        p[end] = s[2] # the last entry of s[1] is zero. we replace it with a d-th root of unity.

        res = track(egtracker, p, 1)
        if is_success(res) && is_finite(res) && is_nonsingular(res)
            q = copy(egtracker.state.solution)
            # we only want nonsingular solutions q
            # an additional check is whether q can be tracked the reverse homotopy from 0 to 1
            # it's enough to go from t = 0 to t = 0.1, since singularities can only be expected at t = 0.
            tracker = egtracker.tracker
            track!(tracker, q, 0.0, 0.1)
            if is_success(status(tracker))
                push!(X, q)
            end
        end
        update_progress_tasks!(progress, i, l_start)
    end

    nothing
end

function threaded_intersection!(X, start, tracker, progress)

    l_start = length(start)
    start_vec = collect(start)

    # Pre-allocate one tracker per thread
    nthr = Threads.nthreads()
    trackers = [deepcopy(tracker) for _ = 1:nthr]

    # Use a lock for thread-safe operations on X and progress
    progress_lock = ReentrantLock()

    next_idx = Threads.Atomic{Int}(1)
    Threads.@sync begin
        for tid = 1:nthr
            let local_egtracker = trackers[tid]
                Threads.@spawn begin
                    while true
                        idx = Threads.atomic_add!(next_idx, 1)
                        if idx > length(start_vec)
                            break
                        end

                        s = start_vec[idx]
                        p = copy(s[1])
                        p[end] = s[2] # the last entry of s[1] is zero. we replace it with a d-th root of unity.

                        res = track(local_egtracker, p, 1)

                        if is_success(res) && is_finite(res) && is_nonsingular(res)
                            q = copy(local_egtracker.state.solution)
                            # we only want nonsingular solutions q
                            # an additional check is whether q can be tracked the reverse homotopy from 0 to 1
                            # it's enough to go from t = 0 to t = 0.1, since singularities can only be expected at t = 0.
                            local_tracker = local_egtracker.tracker
                            track!(local_tracker, q, 0.0, 0.1)
                            if is_success(status(local_tracker))
                                lock(progress_lock) do
                                    push!(X, q)
                                end
                            end
                        end

                        lock(progress_lock) do
                            update_progress_tasks!(progress, idx, l_start)
                        end
                    end
                end
            end
        end
    end

    nothing
end

function set_up_u_homotopy(H, u, W, X, f, g, vars)

    P, Q = get_num_den(g)
    d = degree(P)
    g0 = (u^d - 1) / Q

    # we start with the linear space L which does not use pose conditions on u, so that u^d=1
    # we end with the linear space K with u=c.
    L = linear_subspace_u(W)
    K = linear_subspace(X)

    F₀ = slice(System([f; g0], variables = vars), L; compile = false)
    G₀ = slice(System([f; g], variables = vars), K; compile = false)
    Hom = StraightLineHomotopy(F₀, G₀; gamma = cis(2 * pi * rand()))

    return Hom, d
end



function is_contained(X::WitnessPoints, Y::WitnessSet, F, cache; kwargs...)
    progress = cache.progress
    update_progress!(progress; is_membership_test = true, is_solving = false, is_monodromy = false)
    # main idea: for every x∈X we take a linear space L with codim(L)=dim(Y) through p and move the points in Y to L. Then, we check if the computed points contain x. If yes, return true, else return false.
    LX = linear_subspace(X)
    LY = linear_subspace(Y)

    if length(points(Y)) == 0 || length(points(X)) == 0
        return falses(length(points(X)))
    end
    
    set_up_linear_spaces!(cache, LX, LY)
    out = is_contained(points(X), Y::WitnessSet, F, cache; kwargs...)

    out
end
function is_contained(V::WitnessPoints, W::WitnessSet, cache; kwargs...)
    is_contained(V, W, system(W), cache; kwargs...)
end

function set_up_linear_spaces!(cache, LX, LY)

    n = ambient_dim(LY)
    cX = codim(LX)
    dY = dim(LY)
    k = codim(LY) - codim(LX)

    # to compute linear spaces through the points in X we first set up
    # a matrix-vector pair of the correct size cY×n
    A = cache.A
    b = cache.b

    # since we have used a subset of the equations for Y also for X,
    # we can reuse them
    AX = extrinsic(LX).A
    bX = extrinsic(LX).b

    # the equation for u = 0
    A[1, n] = 1.0
    # new equations
    b[1] = bX[1]
    for i = 2:(k+1)
        for j = 1:n
            A[i, j] = randn(ComplexF64)
        end
    end
    # equations from X, the overlap in linear equations is in the *last* mx-1 equations
    for i = 2:cX
        ℓ = k + i
        for j = 1:n
            A[ℓ, j] = AX[i, j]
        end
        b[ℓ] = bX[i]
    end
end


"""
    DecomposeProgress

"""

Base.@kwdef mutable struct DecomposeProgress
    progress_meter::PM.ProgressUnknown
    current_dim::Int = 0
    n_witness_sets::Int
    n_working_on::Int = 0
    degrees::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    npts::Int = 0
    step::Int = 0
    is_monodromy::Bool = true
    is_finished::Bool = false
end

function update_progress!(progress::DecomposeProgress; is_monodromy = nothing)
    isnothing(is_monodromy) ? nothing : progress.is_monodromy = is_monodromy

    if !isnothing(is_monodromy)
        PM.update!(
            progress.progress_meter,
            progress.step,
            showvalues = showstatus(progress);
            force = true,
        )
    end
end
update_progress!(progress::Nothing, D::Vector{WitnessSet}) = nothing
function update_progress!(progress::DecomposeProgress, W::WitnessSet)
    P = progress.degrees
    c = dim(W)
    if haskey(P, c)
        push!(P[dim(W)], ModelKit.degree(W))
    else
        P[dim(W)] = [ModelKit.degree(W)]
    end

    PM.update!(progress.progress_meter, progress.step, showvalues = showstatus(progress))
end
update_progress_dim!(progress::DecomposeProgress, d::Int) = (progress.current_dim = d)
update_progress_dim!(progress::Nothing, d::Int) = nothing
update_progress_step!(progress::Nothing) = nothing
function update_progress_step!(progress::DecomposeProgress)
    progress.step += 1
end
update_progress_n!(progress::Nothing) = nothing
function update_progress_n!(progress::DecomposeProgress)
    if progress.n_working_on < progress.n_witness_sets
        progress.n_working_on += 1
    end

    PM.update!(progress.progress_meter, progress.step, showvalues = showstatus(progress))
end
update_progress_npts!(progress::Nothing, ℓ::Int) = nothing
function update_progress_npts!(progress::DecomposeProgress, ℓ::Int)
    progress.npts = ℓ

    PM.update!(progress.progress_meter, progress.step, showvalues = showstatus(progress))
end
function finish_progress!(progress::DecomposeProgress)
    progress.is_monodromy = false
    progress.is_finished = true
    progress.current_dim = -1
    PM.finish!(progress.progress_meter, showvalues = showstatus(progress))
end
function showstatus(progress::DecomposeProgress)
    ℓ = progress.npts
    if ℓ > 0
        text = [(
            "Decomposing witness set",
            "$(progress.n_working_on)/$(progress.n_witness_sets) ($ℓ points left to classify)",
        )]
    else
        text = [(
            "Decomposing witness set",
            "$(progress.n_working_on)/$(progress.n_witness_sets)",
        )]
    end
    if !progress.is_finished
        if progress.is_monodromy
            push!(text, ("Status", "running monodromy"))
        else
            push!(text, ("Status", "partitioning points"))
        end
    else
        push!(text, ("Status", "done"))
    end
    degs = progress.degrees
    dims = collect(keys(degs))
    if progress.current_dim ≥ 0 && !in(progress.current_dim, dims)
        push!(dims, progress.current_dim)
    end
    sort!(dims; rev = true)
    for d in dims
        deg = get(degs, d, nothing)
        if isnothing(deg)
            if d == progress.current_dim
                deg_string = "..."
            else
                deg_string = ""
            end
        else
            deg_string = join(map(string, deg), ',')
            if d == progress.current_dim
                deg_string = string(deg_string, ",...")
            end
        end
        push!(text, ("Degrees of components of dim. $d", deg_string))
    end

    text
end


"""
    decompose_with_monodromy!(
        W,
        show_monodromy_progress,
        options,
        max_iters,
        threading,
        progress,
        seed) 

The core function for decomposing a witness set into irreducible components.
"""
function decompose_with_monodromy!(
    W,
    show_monodromy_progress,
    options,
    max_iters,
    warning,
    progress,
    seed;
    threading::Bool = Threads.nthreads() > 1,
)
    update_progress_dim!(progress, dim(W))

    P = points(W)
    L = linear_subspace(W)
    G = system(W)
    n = ambient_dim(L)

    decomposition = Vector{WitnessSet}()

    if dim(L) < n
        update_progress!(progress; is_monodromy = true)

        MS = MonodromySolver(G, L; compile = false, options = options)
        res = monodromy_solve(
            MS,
            P,
            L,
            seed;
            threading = threading,
            show_progress = show_monodromy_progress,
        )
        update_progress!(progress; is_monodromy = false)
        update_progress_npts!(progress, nsolutions(res))

        if warning && (trace(res) > options.trace_test_tol)
            @warn "Trying to decompose non-complete set of witness points for codimension $(dim(L)) (trace test failed). Will try to compute the missing points. The output will contain all components, for which the trace test was successfull."
        end

        iter = 0
        non_complete_points = solutions(res)
        d = length(non_complete_points) # the total degree (i.e., number of points on all irreducible components)
        non_complete_orbits = Vector{Set{Int}}()

        while !isempty(non_complete_points)
            iter += 1
            if iter > max_iters
                break
            end

            if iter > 1
                # for safety an additional monodromy when iter > 1
                update_progress!(progress; is_monodromy = true)
                res = monodromy_solve(
                    MS,
                    non_complete_points,
                    L,
                    seed;
                    threading = threading,
                    show_progress = show_monodromy_progress,
                )
                update_progress!(progress; is_monodromy = false)
            end

            d += nresults(res) - length(non_complete_points) # update total degree in case we found new points.
            non_complete_points = solutions(res)
            ℓ = length(non_complete_points) # once for a later check
            k = length(non_complete_points) # and once for counting progress
            update_progress_npts!(progress, k)


            # Get orbits from monodromy result            
            orbits = get_orbits_from_monodromy_permutations(
                res;
                initial_orbits = non_complete_orbits,
            )

            complete_orbits = Vector{Set{Int}}()

            for orbit in orbits
                update_progress!(progress; is_monodromy = true)

                P_orbit = non_complete_points[collect(orbit)]
                res_orbit = monodromy_solve(
                    MS,
                    P_orbit,
                    L,
                    seed;
                    threading = threading,
                    show_progress = show_monodromy_progress,
                )
                update_progress!(progress; is_monodromy = false)

                if trace(res_orbit) < options.trace_test_tol

                    # We do not want to add orbits of length 1 in the beginning. Even if they are on an irreducible component of degree > 1, they tend to have small trace.
                    if length(orbit) > 1 || iter ≥ 5
                        W_new = WitnessSet(G, L, P_orbit; is_irreducible = true)

                        push!(decomposition, W_new)
                        push!(complete_orbits, orbit)
                        k = k - length(orbit)

                        update_progress!(progress, W_new)
                        update_progress_npts!(progress, k)
                    end
                end
            end

            # Check if we are done
            if sum(degree, decomposition; init = 0) == d
                update_progress_step!(progress)
                break
            end

            # We now have a bunch of partial orbits.
            # For all points in the partial orbits we want to redo a monodromy loop
            # But we also want to reuse the existing orbit information
            # This is a little bit annoyinf since the permutations in our updated monodromy loop will use another set of indices
            # (namely the indices of the points in the partial orbits)
            # Example: We get the orbits (1, 4, 2), (3, 5), (6)
            # Orbit (1, 4, 2) is complete. Non-complete is orbit (3, 5, 6)
            # Then we only do a monodromy loop for the points 3, 5 and 6
            # But the permutation returned from this monodromy loop will be (3, 2) and (1)
            # If we map these indices to the indices of the points in the original set of witness points
            # (i.e. 3 -> 6, 2 -> 5, 1 -> 3)
            # Then we can observe the orbit (3, 5, 6)
            # Here in our implementation we do not actually map the indices to the original set but rather forward, i.e.,
            # (6 -> 3, 5 -> 2, 3 -> 1)

            complete_orbit_indices = Vector{Int}()
            for orbit in complete_orbits
                append!(complete_orbit_indices, collect(orbit))
            end
            sort!(complete_orbit_indices)

            non_complete_points = non_complete_points[setdiff(
                1:length(non_complete_points),
                complete_orbit_indices,
            )]


            orbit_indices_mapping = Dict{Int,Int}()
            # We need to map the orbit indices to the new indices of the non_complete_points
            i = 1
            delta = 0
            for j in complete_orbit_indices
                while i < j
                    orbit_indices_mapping[i] = i - delta
                    i += 1
                end
                # i == j
                delta += 1
                i += 1
            end
            while i <= ℓ
                orbit_indices_mapping[i] = i - delta
                i += 1
            end

            shift_orbit(orbit) = Set(orbit_indices_mapping[i] for i in orbit)
            # 1. shift existing non complete orbits to new mapping
            non_complete_orbits =
                shift_orbit.(
                    merge_sets(
                        [
                            # Remove all orbits that are contained in a complete orbit
                            filter(non_complete_orbits) do o
                                all(co -> isdisjoint(co, o), complete_orbits)
                            end
                            setdiff(orbits, complete_orbits)
                        ],
                    ),
                )
        end
        update_progress_step!(progress)
        update_progress_npts!(progress, 0)
    else
        for p in P
            W_new = WitnessSet(G, L, [p]; is_irreducible = true)
            update_progress!(progress, W_new)
            push!(decomposition, W_new)
        end
        update_progress_step!(progress)
    end

    update_progress!(progress)

    decomposition
end



# Helper function to merge two sets
function merge_sets_find(parent, i)
    if parent[i] == i
        return i
    end
    return merge_sets_find(parent, parent[i])
end
function merge_sets_union(parent, rank, x, y)
    xroot = merge_sets_find(parent, x)
    yroot = merge_sets_find(parent, y)
    if rank[xroot] < rank[yroot]
        parent[xroot] = yroot
    elseif rank[xroot] > rank[yroot]
        parent[yroot] = xroot
    else
        parent[yroot] = xroot
        rank[xroot] += 1
    end
end

# Merge the sets
function merge_sets(sets)
    # Initialize the union-find data structure
    parent = Dict{Int,Int}()
    rank = Dict{Int,Int}()
    for s in sets
        for elem in s
            parent[elem] = elem
            rank[elem] = 0
        end
    end
    # Merge the sets
    for s in sets
        for elem in s
            merge_sets_union(parent, rank, elem, first(s))
        end
    end
    # Build the result
    result = Dict{Int,Set{Int}}()
    for s in sets
        for elem in s
            root = merge_sets_find(parent, elem)
            if !haskey(result, root)
                result[root] = Set{Int}()
            end
            push!(result[root], elem)
        end
    end
    return collect(values(result))
end


# Returns a vector of sets of indices of points in P that are in the same orbit
function get_orbits_from_monodromy_permutations(
    monodromy_result;
    initial_orbits = Vector{Set{Int}}(),
)
    P = permutations(monodromy_result)
    orbits = copy(initial_orbits)

    for col in eachcol(P)
        for (i, j) in enumerate(col)
            # Ignore 0s in orbit since they indicate an error
            if j == 0
                continue
            end
            # Find all orbits where i or j is contained
            matching_orbits = Int[]
            for (k, orbit) in enumerate(orbits)
                if i in orbit || j in orbit
                    push!(matching_orbits, k)
                end
            end


            if (length(matching_orbits) == 1)
                # Add i and j to the orbit
                push!(orbits[matching_orbits[1]], i, j)
            else
                # Merge all matching orbits
                merged_orbit = Set([i, j])
                for k in reverse(matching_orbits)
                    union!(merged_orbit, orbits[k])
                    deleteat!(orbits, k)
                end
                push!(orbits, merged_orbit)
            end
        end
    end
    orbits
end

function decompose_with_monodromy_options(M::MonodromyOptions)

    MonodromyOptions(;
        permutations = true,
        trace_test = true,
        single_loop_per_start_solution = true,
        check_startsolutions = M.check_startsolutions,
        group_actions = M.group_actions,
        loop_finished_callback = M.loop_finished_callback,
        parameter_sampler = M.parameter_sampler,
        equivalence_classes = M.equivalence_classes,
        trace_test_tol = M.trace_test_tol,
        target_solutions_count = M.target_solutions_count,
        timeout = M.timeout,
        min_solutions = M.min_solutions,
        max_loops_no_progress = M.max_loops_no_progress,
        reuse_loops = M.reuse_loops,
        distance = M.distance,
        unique_points_atol = M.unique_points_atol,
        unique_points_rtol = M.unique_points_rtol,
    )
end


"""
    decompose(Ws::Vector{WitnessPoints}; options...) 

This function decomposes a [`WitnessSet`](@ref) or a vector of [`WitnessSet`](@ref) into irreducible components.

### Options
* `show_progress = true`: indicate whether a progress bar should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) should be displayed.
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_iters = 50`: maximal number of iterations for the decomposition step.
* `warning = true`: if `true` prints a warning when the [`trace_test`](@ref) fails. 
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)
* `seed`: choose the random seed.


### Example
The following example decomposes the witness set for a union of two circles.

```julia-repl
julia> @var x y
julia> f = (x^2 + y^2 - 1) * ((x-1)^2 + (y-1)^2 - 1)
julia> W = regeneration([f])
julia> dec = decompose(W)
2-element Vector{WitnessSet}:
 Witness set for dimension 1 of degree 2
 Witness set for dimension 1 of degree 2
```

The decomposition can be stored in a [`NumericalIrreducibleDecomposition`](@ref) data structure.
```julia-repl
julia> N = NumericalIrreducibleDecomposition(dec)
Numerical irreducible decomposition with 2 components
=====================================================
• 2 component(s) of dimension 1.

 degree table of components:
╭───────────┬───────────────────────╮
│ dimension │ degrees of components │
├───────────┼───────────────────────┤
│     1     │        (2, 2)         │
╰───────────┴───────────────────────╯
```

"""
function decompose(
    Ws::Union{Vector{WP}};
    show_progress::Bool = true,
    show_monodromy_progress::Bool = false,
    monodromy_options::MonodromyOptions = MonodromyOptions(;
        trace_test_tol = 1e-10,
        parameter_sampler = weighted_normal,
    ),
    max_iters::Int = 50,
    warning::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    seed = nothing,
) where {WP<:WitnessSet}

    if isnothing(seed)
        seed = rand(UInt32)
    end
    Random.seed!(seed)

    sort!(Ws; by = dim, rev = true)
    options = decompose_with_monodromy_options(monodromy_options)
    out = Vector{WitnessSet}()

    if isempty(Ws)
        return out
    end

    c = length(Ws)
    n = ambient_dim(linear_subspace(Ws[1]))

    if show_progress
        progress = DecomposeProgress(
            progress_meter = PM.ProgressUnknown(
                dt = 0.1,
                desc = "Decomposing $c witness sets",
                enabled = true,
                spinner = true,
            ),
            n_witness_sets = c,
        )
    else
        progress = nothing
    end

    update_progress_n!(progress)
    update_progress_step!(progress)

    for W in Ws
        if ModelKit.degree(W) > 0

            dec = decompose_with_monodromy!(
                W,
                show_monodromy_progress,
                options,
                max_iters,
                warning,
                progress,
                seed;
                threading = threading,
            )
            if !isnothing(dec)
                append!(out, dec)
            end

            update_progress_n!(progress)
        end

    end
    finish_progress!(progress)

    out
end
decompose(W::WitnessSet; kwargs...) = decompose([W]; kwargs...)



"""
    NumericalIrreducibleDecomposition

Store the witness sets in a common data structure.
"""
struct NumericalIrreducibleDecomposition
    Witness_Sets::Dict
    seed::Any
end
NumericalIrreducibleDecomposition(Ws::Vector{WitnessSet}) =
    NumericalIrreducibleDecomposition(Ws, nothing)
function NumericalIrreducibleDecomposition(Ws::Vector{WitnessSet}, seed)
    D = Dict{Int,Vector{WitnessSet}}()
    for W in Ws
        k = dim(W)
        if !haskey(D, k)
            D[k] = [W]
        else
            push!(D[k], W)
        end
    end
    NumericalIrreducibleDecomposition(D, seed)
end

"""
    witness_sets(N::NumericalIrreducibleDecomposition;
        dims::Union{Vector{Int},Nothing} = nothing)

Returns the witness sets in `N`. 
`dims` specifies the dimensions that should be considered.
"""
function witness_sets(
    N::NumericalIrreducibleDecomposition;
    dims::Union{Vector{Int},Nothing} = nothing,
)
    D = N.Witness_Sets
    if isnothing(dims)
        out = D
    else
        out = Dict{Int,Vector{WitnessSet}}()
        for k in dims
            if haskey(D, k)
                out[k] = D[k]
            end
        end
    end

    out
end
witness_sets(N::NumericalIrreducibleDecomposition, dim::Int) = witness_sets(N, [dim])
seed(N::NumericalIrreducibleDecomposition) = N.seed

"""
    ncomponents(N::NumericalIrreducibleDecomposition;
        dims::Union{Vector{Int},Nothing} = nothing)

Returns the total number of components in `N`. 
`dims` specifies the dimensions that should be considered.
"""
function ncomponents(
    N::NumericalIrreducibleDecomposition;
    dims::Union{Vector{Int},Nothing} = nothing,
)
    D = witness_sets(N)
    if isempty(D)
        return 0
    elseif isnothing(dims)
        out = sum(length(last(Ws)) for Ws in D)
    else
        out = 0
        for d in dims
            if haskey(D, d)
                out += length(D[d])
            end
        end
    end

    out
end
ncomponents(N::NumericalIrreducibleDecomposition, dim::Int) = ncomponents(N; dims = [dim])
n_components(N; dims = nothing) = ncomponents(N; dims = dims)
n_components(N, dim) = ncomponents(N; dims = [dim])

"""

    degrees(N::NumericalIrreducibleDecomposition;
        dims::Union{Vector{Int},Nothing} = nothing)

Returns the degrees of the components in `N`.
`dims` specifies the dimensions that should be considered.

"""
function ModelKit.degrees(
    N::NumericalIrreducibleDecomposition;
    dims::Union{Vector{Int},Nothing} = nothing,
)
    D = N.Witness_Sets
    out = Dict{Int,Vector{Int}}()
    if isnothing(dims)
        for key in keys(D)
            out[key] = degree.(D[key])
        end
    else
        out = Dict{Int,Vector{Int}}()
        if !isempty(dims)
            for k in dims
                if haskey(D, k)
                    out[k] = degree.(D[k])
                end
            end
        end
    end

    out
end
ModelKit.degrees(N::NumericalIrreducibleDecomposition, dim::Int) = degrees(N, [dim])

function max_dim(N::NumericalIrreducibleDecomposition)
    D = witness_sets(N)
    k = keys(D)
    if !isempty(k)
        maximum(k)
    else
        -1
    end
end

function Base.show(io::IO, N::NumericalIrreducibleDecomposition)
    D = witness_sets(N)
    if !isempty(D)
        total = sum(length(last(Ws)) for Ws in D)
    else
        total = 0
    end
    if total == 1
        header = "Numerical irreducible decomposition with 1 component"
    else
        header = "Numerical irreducible decomposition with $total components"
    end
    println(io, header)
    println(io, "="^(length(header)))
    mdim = max_dim(N)
    if mdim >= 0
        for d = max_dim(N):-1:0
            if haskey(D, d)
                ℓ = length(D[d])
                if ℓ > 0
                    println(io, "• $ℓ component(s) of dimension $d.")
                end
            end
        end

        println(io, "\n degree table of components:")
        degree_table(io, N)
    end
end
function degree_table(io, N::NumericalIrreducibleDecomposition)
    D = witness_sets(N)
    k = collect(keys(D))
    sort!(k, rev = true)
    n = length(k)

    headers = ["dimension", "degrees of components"]
    data = Matrix{Union{Int,String}}(undef, n, 2)

    for (i, key) in enumerate(k)
        data[i, 1] = key
        components = [ModelKit.degree(W) for W in D[key]]
        sort!(components, rev = true)
        if length(components) == 1
            data[i, 2] = first(components)
        elseif length(components) <= 10
            s = string(components)
            data[i, 2] = string("(", s[2:(end-1)], ")")
        else
            s = string(components[1:10])
            data[i, 2] = string("(", s[2:(end-1)], ", ...)")
        end
    end

    PrettyTables.pretty_table(
        io,
        data;
        column_labels = headers,
        backend = :text,
        table_format = PrettyTables.TextTableFormat(
            borders = PrettyTables.text_table_borders__unicode_rounded,
        ),
        alignment = :c,
        style = PrettyTables.TextTableStyle(
            column_label = PrettyTables.Crayon(bold = false),
            table_border = PrettyTables.Crayon(faint = true),
        ),
    )
end

"""
    numerical_irreducible_decomposition(F::System; options...)

Computes the numerical irreducible of the variety defined by ``F=0``. 

### Options

* `show_progress = true`: indicate whether a progress bar should be displayed.
* `sorted = true`: the polynomials in F will be sorted by degree in decreasing order. 
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `monodromy_options_for_regeneration`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref) in [`regeneration`](@ref).
* `monodromy_options_for_decompose`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref) in [`decompose`](@ref).
* `show_progresss = false`: if `true`, sets `show_monodromy_for_regeneration_progress` and `show_monodromy_for_decompose_progress` to `true`.
* `show_monodromy_for_regeneration_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) in [`regeneration`](@ref) should be displayed.
* `show_monodromy_for_decompose_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) in [`decompose`](@ref) should be displayed.
* `max_iters = 50`: maximal number of iterations for the decomposition step.
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).` This option is used for [`regeneration`](@ref).
* `warning = true`: if `true` prints a warning when the [`trace_test`](@ref) fails. 
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)
* `seed`: choose the random seed.

### Example
The following example computes witness sets for a union of one 2-dimensional component of degree 2, two 1-dimensional components each of degree 4 and 8 points.
```julia-repl
julia> @var x, y, z
julia> p = (x * y - x^2) + 1 - z
julia> q = x^4 + x^2 - y - 1
julia> F = [p * q * (x - 3) * (x - 5);
            p * q * (y - 3) * (y - 5);
            p * (z - 3) * (z - 5)]

julia> N = numerical_irreducible_decomposition(F)
Numerical irreducible decomposition with 11 components
======================================================
• 1 component(s) of dimension 2.
• 2 component(s) of dimension 1.
• 8 component(s) of dimension 0.

 degree table of components:
╭───────────┬──────────────────────────╮
│ dimension │  degrees of components   │
├───────────┼──────────────────────────┤
│     2     │            2             │
│     1     │          (4, 4)          │
│     0     │ (1, 1, 1, 1, 1, 1, 1, 1) │
╰───────────┴──────────────────────────╯
```

`numerical_irreducible_decomposition` also works for systems of rational functions:
```julia-repl
julia> @var x y z
julia> g = [x^2 + y^2 - z; x/(y-1) + y + z - 1]
julia> N = numerical_irreducible_decomposition(g)
Numerical irreducible decomposition with 1 component
====================================================
• 1 component(s) of dimension 1.

 degree table of components:
╭───────────┬───────────────────────╮
│ dimension │ degrees of components │
├───────────┼───────────────────────┤
│     1     │           4           │
╰───────────┴───────────────────────╯
```
"""
function numerical_irreducible_decomposition(
    F::System;
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_accuracy = 1e-10,
    ),
    monodromy_options_for_regeneration = MonodromyOptions(;
        trace_test = true,
        parameter_sampler = weighted_normal,
    ),
    monodromy_options_for_decompose::MonodromyOptions = MonodromyOptions(;
        trace_test_tol = 1e-10,
    ),
    show_monodromy_progress::Bool = false,
    show_monodromy_for_regeneration_progress::Bool = false,
    show_monodromy_for_decompose_progress::Bool = false,
    max_iters::Int = 50,
    sorted::Bool = true,
    max_codim::Union{Int,Nothing} = nothing,
    warning::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    seed = nothing,
    atol = 1e-14,
    rtol = sqrt(eps()),
    kwargs...,
)
    if show_monodromy_progress
        show_monodromy_for_regeneration_progress = true
        show_monodromy_for_decompose_progress = true
    end

    !isnothing(seed) && Random.seed!(seed)

    Ws = regeneration(
        F;
        sorted = sorted,
        max_codim = max_codim,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
        monodromy_options = monodromy_options_for_regeneration,
        show_monodromy_progress = show_monodromy_for_regeneration_progress,
        threading = threading,
        seed = nothing,
        atol = atol,
        rtol = rtol,
        kwargs...,
    )
    if isnothing(Ws)
        return nothing
    end
    dec = decompose(
        Ws;
        monodromy_options = monodromy_options_for_decompose,
        max_iters = max_iters,
        show_monodromy_progress = show_monodromy_for_decompose_progress,
        threading = threading,
        warning = warning,
        seed = seed,
        kwargs...,
    )

    NumericalIrreducibleDecomposition(dec, seed)
end
"""
    nid(F::System; options...)

Calls [`numerical_irreducible_decomposition`](@ref).

### Example
```julia-repl
julia> @var x, y
julia> f = [x^2 + y^2 - 1]
julia> N = nid(f)
Numerical irreducible decomposition with 1 component
====================================================
• 1 component(s) of dimension 1.

 degree table of components:
╭───────────┬───────────────────────╮
│ dimension │ degrees of components │
├───────────┼───────────────────────┤
│     1     │           2           │
╰───────────┴───────────────────────╯
```
"""
nid(F::System; kwargs...) = numerical_irreducible_decomposition(F; kwargs...)
nid(F::Vector{Expression}; kwargs...) = numerical_irreducible_decomposition(F; kwargs...)
nid(F::Expression; kwargs...) = numerical_irreducible_decomposition([F]; kwargs...)
numerical_irreducible_decomposition(F::Vector{Expression}; kwargs...) =
    numerical_irreducible_decomposition(System(F); kwargs...)


### intersecting witness sets with hypersurfaces 
"""
    IntersectProgress

"""
Base.@kwdef mutable struct IntersectProgress
    progress_meter::PM.ProgressUnknown
    is_solving::Bool = false
    is_membership_test::Bool = false
    is_monodromy::Bool = false
    is_finished::Bool = false
    current_task::Int = 0
    ntasks::Int = 0
end
IntersectProgress(progress_meter::PM.ProgressUnknown) = IntersectProgress(progress_meter = progress_meter)
function update_progress!(
    progress::IntersectProgress;
    is_solving::Union{Nothing,Bool} = nothing,
    is_membership_test::Union{Nothing,Bool} = nothing,
    is_monodromy::Union{Nothing,Bool} = nothing,
)
    isnothing(is_solving) ? nothing : progress.is_solving = is_solving
    isnothing(is_membership_test) ? nothing :
    progress.is_membership_test = is_membership_test
    isnothing(is_monodromy) ? nothing : progress.is_monodromy = is_monodromy
    PM.update!(
        progress.progress_meter,
        showvalues = showvalues(progress);
        force = true,
    )
end
update_progress!(progress::IntersectProgress, W) = nothing
function update_progress_tasks!(progress::IntersectProgress, i::Int, m::Int)
    progress.current_task = i
    progress.ntasks = m
    PM.update!(
        progress.progress_meter,
        showvalues = showvalues(progress),
    )
end
function finish_progress!(progress::IntersectProgress)
    progress.is_solving = false
    progress.is_membership_test = false
    progress.is_monodromy = false
    progress.is_finished = true
    PM.finish!(progress.progress_meter, showvalues = showvalues(progress))
end
function showvalues(progress::IntersectProgress)

    if !progress.is_finished
        text = [("Status", "")]
        if progress.is_solving
            push!(text, ("Track paths", "$(progress.current_task)/$(progress.ntasks)"))
        elseif progress.is_monodromy
            push!(
                text,
                (
                    "Track paths",
                    "$(progress.ntasks)/$(progress.ntasks) (fill up points...)",
                ),
            )
        elseif progress.is_membership_test
            push!(
                text,
                (
                    "Membership test",
                    "$(progress.current_task)/$(progress.ntasks) points checked",
                ),
            )
        end
    else
        text = [("Status", "done")]
    end

    text
end

"""
    IntersectCache

A cache for [`intersect`](@ref).
"""
mutable struct IntersectCache{Sys<:AbstractSystem}
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    x0::Vector{ComplexF64}
    y0::Vector{ComplexF64}
    y::Vector{ComplexF64}

    Fᵢ::Sys
    u::Variable

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{IntersectProgress,Nothing}
end

function IntersectCache(u, F, W₁, W₂, Hᵤ, EO, TO, progress)
    m, N = size(F)
    A = zeros(ComplexF64, N - 1, N)
    b = zeros(ComplexF64, N - 1)
    x0 = zeros(ComplexF64, N)
    y0 = zeros(ComplexF64, m)
    y = zeros(ComplexF64, m)

    IntersectCache(A, b, x0, y0, y, F, u, EO, TO, progress)
end


"""
    intersect

"""
function Base.intersect(W::WitnessSet, H::WitnessSet; 
                    show_progress::Bool = true,
                    tracker_options = TrackerOptions(),
                    endgame_options = EndgameOptions(;
                            max_endgame_steps = 100,
                            max_endgame_extended_steps = 100,
                            sing_cond = 1e12,
                        ),
                    monodromy_options::MonodromyOptions = MonodromyOptions(;
                            trace_test = true,
                            parameter_sampler = weighted_normal,
                        ),
                    show_monodromy_progress::Bool = false,
                    threading = Threads.nthreads() > 1,
                    kwargs...)
    @assert size(system(H), 1) == 1 "The second argument must be defined by a single polynomial."
    @assert size(system(W), 2) == size(system(H), 2) "Witness sets must be in the same ambient space."

    # progress bar
    if show_progress
        progress = IntersectProgress(
            PM.ProgressUnknown(
                dt = 0.001,
                desc = "Intersecting...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end
        

    # transform W and H so that they use the additional variable u
    n = ambient_dim(W.L)
    @unique_var u
    @unique_var vars[1:n]
    W₁, W₂, Hᵤ, FW, fW, fH, vars_u = prepare_for_u_homotopy(H, W, vars, u)

    # system for W ∩ H
    G_sys = System([fW; fH], variables = vars_u)
    G_u = fixed(G_sys; compile = false)

    # cache
    cache = IntersectCache(u, FW, W₁, W₂, Hᵤ, endgame_options, tracker_options, progress)

    
    # intersect
    intersect_with_hypersurface!(W₁, Hᵤ, W₂, cache; threading = threading, kwargs...)
    update_Fᵢ!(cache, G_u)
    fill_up!([W₁; W₂],
            monodromy_options,
            cache,
            show_monodromy_progress,
            threading,
        )

    # return data 
    G = fixed(System([fW; fH], variables = vars); compile = false)
    P1, L1 = u_transform(W₁)
    out = [WitnessSet(G, L1, P1)]
    if !isnothing(W₂)
        P2, L2 = u_transform(W₂)
        push!(out, WitnessSet(deepcopy(G), L2, P2))
        
    end
    filter!(X -> degree(X) > 0, out)

    finish_progress!(progress)
    if length(out) == 1
        return first(out)
    else
        return out
    end
end

function prepare_for_u_homotopy(H, W, vars, u)
    # prepare polynomials 
    FH0 = System(system(H)) |> deepcopy
    FW0 = System(system(W)) |> deepcopy
    fH = subs(expressions(FH0), variables(FH0) => vars)
    fW = subs(expressions(FW0), variables(FW0) => vars)
    vars_u = [vars; u]
    FH = System(fH, variables = vars_u)
    FW = System(fW, variables = vars_u)

    # prepare flags
    LW = linear_subspace(W)
    LH = linear_subspace(H)
    flagW = get_flag(1:2, LW) # returns L, Lᵤ
    flagH = get_flag(1:1, LH)
    cW, cH = extrinsic((flagW[1][1])).b[1], extrinsic((flagH[1][1])).b[1] # the first entry of b is c in "u=c"

    # prepare WitnessPoints
    # in u-coordiantes the points in W have to be transformed as x -> [x; c]
    W₁ = WitnessPoints(flagW[1][1], flagW[1][2], map(x -> [x; cW], solutions(W)))
    if dim(W) == 0 
        W₂ = nothing
    else
        W₂ = WitnessPoints(flagW[2][1], flagW[2][2], Vector{Vector{ComplexF64}}())
    end
    Hᵤ = WitnessSet(fixed(FH; compile = false), flagH[1][1],  map(x -> [x; cH], solutions(H)))

    W₁, W₂, Hᵤ, fixed(FW; compile = false), fW, fH, vars_u
end
