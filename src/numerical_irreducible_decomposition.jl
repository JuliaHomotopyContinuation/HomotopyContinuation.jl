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
linear_subspace(W::WitnessPoints) = W.L
linear_subspace_u(W::WitnessPoints) = W.Lᵤ
points(W::WitnessPoints) = W.R
codim(W::WitnessPoints) = dim(linear_subspace(W))
dim(W::WitnessPoints) = codim(linear_subspace(W))
ModelKit.degree(W::WitnessPoints) = length(points(W))
push!(W::WitnessPoints, R::Vector{ComplexF64}) = push!(W.R, R)

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
    path_label::String = "Track paths"
    current_path::Int = 0
    npaths::Int = 0
    current_task::Int = 0
    ntasks::Int = 0
    current_hypersurface::Int = 1
    nhypersurfaces::Int
    monodromy_codim::Int = 0
    monodromy_solutions::Int = 0
    monodromy_tracked_loops::Int = 0
    monodromy_generated_loops::Int = 0
    monodromy_no_change::Int = 0
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
    )
end
function update_progress!(
    progress::WitnessSetsProgress,
    stats::MonodromyStatistics;
    queued::Int,
    solutions::Int,
    finish::Bool = false,
)
    progress.monodromy_solutions = solutions
    progress.monodromy_tracked_loops = stats.tracked_loops[]
    progress.monodromy_generated_loops = stats.generated_loops[]
    progress.monodromy_no_change = loops_no_change(stats, solutions)

    if finish
        PM.update!(
            progress.progress_meter,
            progress.current_hypersurface,
            showvalues = showvalues(progress),
        )
    elseif time() > progress.progress_meter.tlast + progress.progress_meter.dt
        PM.update!(
            progress.progress_meter,
            progress.current_hypersurface,
            showvalues = showvalues(progress),
        )
    end

    nothing
end
update_progress_monodromy_witness!(progress::Nothing, W) = nothing
function update_progress_monodromy_witness!(progress::WitnessSetsProgress, W)
    progress.monodromy_codim = codim(W)
    progress.monodromy_solutions = degree(W)
    progress.monodromy_tracked_loops = 0
    progress.monodromy_generated_loops = 0
    progress.monodromy_no_change = 0
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress_hypersurface!(progress::Nothing, i::Int) = nothing
function update_progress_hypersurface!(progress::WitnessSetsProgress, i::Int)
    progress.current_hypersurface = i
    progress.current_path = 0
    progress.npaths = 0
    progress.path_label = "Track paths"
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
update_progress_paths!(progress::Nothing, i::Int, m::Int) = nothing
function update_progress_paths!(progress::WitnessSetsProgress, i::Int, m::Int)
    progress.current_path = i
    progress.npaths = m
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress_path_label!(progress::Nothing, label::String) = nothing
function update_progress_path_label!(progress::WitnessSetsProgress, label::String)
    progress.path_label = label
    nothing
end
start_progress_paths!(progress::Nothing, label::String, m::Int) = nothing
function start_progress_paths!(progress::WitnessSetsProgress, label::String, m::Int)
    progress.path_label = label
    progress.current_path = 0
    progress.npaths = m
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
            if progress.is_solving && progress.npaths > 0
                push!(
                    text,
                    (progress.path_label, "$(progress.current_path) / $(progress.npaths)"),
                )
            elseif progress.is_monodromy
                push!(
                    text,
                    (
                        "Fill up points",
                        "$(progress.monodromy_tracked_loops) loops tracked ($(progress.monodromy_generated_loops) generated, $(progress.monodromy_no_change) no change)",
                    ),
                )
            elseif progress.is_membership_test
                push!(
                    text,
                    ("Membership test", "$(progress.current_task) / $(progress.ntasks)"),
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
        if progress.is_monodromy && c == progress.monodromy_codim
            deg = progress.monodromy_solutions
        end
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
    # array caches
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    x0::Vector{ComplexF64}
    y0::Vector{ComplexF64}
    y::Vector{ComplexF64}

    # carry polynomials along
    fᵢ::Vector{Expression}
    Fᵢ::Sys
    h::Expression
    u::Variable
    i::Int
    codim::Int
    projective::Bool
    ℓ::Union{Nothing,Expression}
    ℓ_coeffs::Union{Nothing,Vector{ComplexF64}}
    max_trials_u_homotopy::Int
    intrinsic::Union{Nothing,Bool}

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{WitnessSetsProgress,Nothing}
end

function RegenerationCache(
    u,
    f,
    F,
    h,
    codim,
    projective,
    ℓ,
    ℓ_coeffs,
    max_trials_u_homotopy,
    intrinsic,
    EO,
    TO,
    progress,
)
    m, N = size(F)
    A = zeros(ComplexF64, N - 1, N)
    b = zeros(ComplexF64, N - 1)
    x0 = zeros(ComplexF64, N)
    y0 = zeros(ComplexF64, m)
    y = zeros(ComplexF64, m)

    RegenerationCache(
        A,
        b,
        x0,
        y0,
        y,
        f,
        F,
        h,
        u,
        0,
        codim,
        projective,
        ℓ,
        ℓ_coeffs,
        max_trials_u_homotopy,
        intrinsic,
        EO,
        TO,
        progress,
    )
end
function update_Fᵢ!(cache, fᵢ, vars)
    Fᵢ = fixed(System(fᵢ, variables = vars), compile = false)
    cache.fᵢ = fᵢ
    cache.Fᵢ = Fᵢ
    m = length(Fᵢ)
    resize!(cache.y0, m)
    fill!(cache.y0, zero(ComplexF64))
    resize!(cache.y, m)
    fill!(cache.y, zero(ComplexF64))

    nothing
end
update_hᵢ!(cache, h) = cache.h = h
update_i!(cache, i) = cache.i = i


"""
    regeneration(F::System; options...) 

This solves ``F=0`` equation-by-equations and returns a [`WitnessSet`](@ref) for every dimension without decomposing them into irreducible components (witness sets that are not decomposed are also called witness supersets).

The implementation is based on the algorithm [u-regeneration](https://arxiv.org/abs/2206.02869) by Duff, Leykin and Rodriguez. 

### Options
* `sorted`: if `true` (default), the polynomials in `F`` will be sorted by degree in increasing order (when `F` is a system of rational functions, we first compute witness sets for each entry of `F` and then sort `F` according to their degrees). If `false, the polynomials in `F` will not be sorted. The third option is `:randomized`, which multiplies `F` by a random matrix.
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `show_progress = true`: indicate whether a progress bar should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) should be displayed. If `false`, minimal info about the monodromy computations are still displayed in the progress bar of `regeneration`.
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_trials_u_homotopy = 5`: maximal number of random subspaces tried until an intermediate u-regeneration intersection step succeeds.
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).`
* `warning = true`: if `true`, prints a warning when an intermediate u-regeneration intersection step fails after all trials.
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
    sorted::Union{Bool,Symbol} = true,
    max_codim::Union{Int,Nothing} = nothing,
    show_progress::Bool = true,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_cond = 1e12,
    ),
    monodromy_options::MonodromyOptions = MonodromyOptions(),
    max_trials_u_homotopy::Int = 5,
    intrinsic::Union{Nothing,Bool} = true, # always intrinsic linear subspace homotopies
    show_monodromy_progress::Bool = false,
    threading::Bool = Threads.nthreads() > 1,
    seed = nothing,
    atol = 1e-14,
    rtol = sqrt(eps()),
    warning::Bool = true,
)

    !isnothing(seed) && Random.seed!(seed)

    # the algorithm is u-regeneration as proposed 
    # by Duff, Leykin and Rodriguez in https://arxiv.org/abs/2206.02869

    # declare variables:
    # u-regeneration adds another variable u to F
    xvars = variables(F)
    projective = is_homogeneous(F)
    vars = copy(xvars)
    @unique_var u
    push!(vars, u)

    n = size(F, 2) # ambient dimension
    c = size(F, 1) # we can have witness sets of codimesion at most min(c,n)
    expected_max_codim = min(c, n - projective)
    if !isnothing(max_codim) && max_codim < expected_max_codim
        # if max_codim is smaller than the expected codimension we must compute witness points for one more codimension, so that we can remove spurious points
        codim = max_codim + 1
    else
        codim = expected_max_codim
    end

    # prepare equations
    f = expressions(F)
    is_poly = all(is_polynomial, f)
    if is_poly && prepare_polynomials!(f, sorted, projective)
        return nothing
    end

    # progress bar
    if show_progress
        progress_meter = PM.ProgressUnknown(
            dt = 0.4,
            desc = "Computing witness sets...",
            enabled = true,
            spinner = true,
        )
        progress_meter.tlast += 0.3
        progress = WitnessSetsProgress(n, c, progress_meter)
    else
        progress = nothing
    end

    # prepare codim witness sets for the output
    # internally we represent a witness superset by WitnessPoints
    # the i-th witness superset out[i] is for codimension i
    out = initialize_witness_sets(codim, n; affine = !projective)

    # we compute witness (super)sets for the hypersurfaces f[1]=0,...,f[c]=0.
    # it is covenient to use the WitnessSet wrapper here, because this also keeps track of the equation
    # as a linear subspace we take the linear subspace for out[1], that sets u=0.
    update_progress!(progress; is_computing_hypersurfaces = true)
    H = initialize_hypersurfaces(
        f,
        vars,
        xvars,
        linear_subspace(out[1]),
        projective;
        threading = threading,
    )
    if any(isnothing, H)
        return nothing
    end

    # # sort expressions by degree of hypersurfaces when we have rational functions 
    if sorted == true && !is_poly
        σ = sortperm(H, by = ModelKit.degree)
        f = f[σ]
        H = H[σ]
    end

    # Initialize a cache
    Fᵢ = fixed(System(f[1:1], variables = vars), compile = false)
    ℓ_coeffs = projective ? randn(ComplexF64, length(xvars)) : nothing
    ℓ = projective ? sum(ℓ_coeffs .* xvars) : nothing
    cache = RegenerationCache(
        u,
        f[1:1],
        Fᵢ,
        f[1],
        codim,
        projective,
        ℓ,
        ℓ_coeffs,
        max_trials_u_homotopy,
        intrinsic,
        endgame_options,
        tracker_options,
        progress,
    )

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
                    update_hᵢ!(cache, f[i])

                    # for all W in out we intersect W with H[i]
                    intersect_all!(
                        out,
                        H,
                        cache;
                        threading = threading,
                        atol = atol,
                        rtol = rtol,
                        warning = warning,
                    )

                    # update Fᵢ
                    update_Fᵢ!(cache, f[1:i], vars)

                    # running one monodromy per witness set to fill up point
                    fill_up!(
                        out,
                        monodromy_options,
                        cache,
                        show_monodromy_progress,
                        threading;
                        atol = atol,
                        rtol = rtol,
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
        ws = map(out) do W
            update_progress!(progress)
            P, L = u_transform(W)
            WitnessSet(fixed(F, compile = false), L, P; projective = projective)
        end
    else
        return ws = Vector{WitnessSet}()
    end


    finish_progress!(progress)

    return ws
end

function prepare_polynomials!(f, sorted, projective)
    if sorted == true
        sort!(f, by = ModelKit.degree)
    elseif sorted == :randomized
        if projective
            @error "Randomization is not available for homogeneous systems."
            return true
        else
            sort!(f, by = ModelKit.degree, rev = true)
            # random triangular system
            c = length(f)
            g = map(1:c) do i
                sum(randn(ComplexF64) * f[j] for j = i:c)
            end
            f .= g
        end
    end
    return false
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
    # type 2 sets the linear equation u=c. If L₀ is linear, keep the augmented
    # subspaces linear as well.
    c = is_linear(L₀) ? zero(ComplexF64) : randn(ComplexF64)
    A = [zeros(1, m) 1.0; A₀ zeros(n - 1)]
    b = [c; b₀]

    map(iter) do i
        # type 1 linear equation for out[i] takes the last i rows of Aᵤ
        idx = i <= n - 1 ? (i:n-1) : (n:n-1)
        Eᵤ = ExtrinsicDescription(Aᵤ[idx, :], bᵤ[idx]; orthonormal = true)
        Lᵤ = LinearSubspace(Eᵤ)
        # type 2 linear equation for out[i] takes the last i-1 rows and the first row of A
        rows = [1; idx .+ 1]
        E = ExtrinsicDescription(A[rows, :], b[rows]; orthonormal = true)
        L = LinearSubspace(E)
        L, Lᵤ
    end
end

function initialize_witness_sets(codim, n; affine::Bool = true)
    # we need an initial slice from which the u-regeneration flag is derived
    dim = affine ? 1 : 2
    L₀ = rand_subspace(n; dim = dim, affine = affine)
    flag = get_flag(1:codim, L₀)

    map(flag) do F
        L, Lᵤ = F
        WitnessPoints(L, Lᵤ, Vector{Vector{ComplexF64}}())
    end
end
function initialize_hypersurfaces(
    f::Vector{Expression},
    vars,
    xvars,
    L,
    projective;
    threading::Bool = Threads.nthreads() > 1,
)
    c = length(f)
    out = Vector{WitnessSet}(undef, c)
    for i = 1:c
        fᵢ = f[i]
        h = fixed(System([fᵢ], variables = vars), compile = false)
        pᵢ, qᵢ = get_num_den(fᵢ) # fᵢ = pᵢ / qᵢ
        if projective
            E = extrinsic(L)
            A = E.A[2:end, 1:(end-1)]
            b = E.b[2:end]
            L_eqs = map(1:size(A, 1)) do i
                sum(A[i, j] * xvars[j] for j = 1:length(xvars)) - b[i]
            end
            res = solve(
                System([pᵢ; L_eqs], variables = xvars);
                start_system = :total_degree,
                compile = false,
                show_progress = false,
                threading = threading,
            )
        else
            res = solve(
                System([pᵢ], variables = vars),
                target_subspace = L;
                start_system = :total_degree,
                compile = false,
                show_progress = false,
                threading = threading,
            )
        end
        # if fᵢ is rational we check which zeros of pᵢ are also zeros of fᵢ
        if qᵢ != 1
            res = solve(
                h,
                projective ? map(x -> [x; 0], solutions(res)) : solutions(res);
                start_subspace = L,
                target_subspace = L,
                show_progress = false,
                threading = threading,
            )
        end
        if isnothing(res)
            return nothing
        end
        if projective && qᵢ == 1
            S = map(x -> [x; 0], solution.(results(res)))
        else
            S = solution.(results(res))
        end
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

    for k = length(out):-1:1  # k = codim(W) for W in Ws
        Wₖ = out[k]
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

function fill_up!(
    out,
    monodromy_options,
    cache,
    show_monodromy_progress,
    threading;
    atol = 1e-14,
    rtol = sqrt(eps()),
)
    progress = cache.progress
    Fᵢ = cache.Fᵢ

    update_progress!(
        progress;
        is_membership_test = false,
        is_solving = false,
        is_monodromy = true,
    )
    for W in out
        if !isnothing(W) && dim(W) > 0 && degree(W) > 0
            update_progress_monodromy_witness!(progress, W)
            L = linear_subspace(W)
            MS = MonodromySolver(
                Fᵢ,
                L;
                options = regeneration_monodromy_options(
                    monodromy_options,
                    W;
                    atol = atol,
                    rtol = rtol,
                ),
            )
            # the witness sets are usually complete already: certify this with the trace before running a full monodromy loop
            if trace_precheck!(MS, W.R, L, threading)
                sols = W.R
                W.R =
                    length(sols) > 1 ? unique_points(sols; atol = atol, rtol = rtol) : sols
            else
                res = monodromy_solve(
                    MS,
                    W.R,
                    L,
                    rand(UInt32);
                    show_progress = show_monodromy_progress,
                    progress = show_monodromy_progress ? nothing : progress,
                    threading = threading,
                    warning = false,
                )
                # check if solutions remain
                if nsolutions(res) == 0
                    W.R = Vector{Vector{ComplexF64}}()
                else
                    sols = solution.(results(res))
                    # Near-singular solutions tracked in monodromy loops can have poor accuracy (~1e-10 error).
                    # If they evade monodromy's internal deduplication tolerance, remove near-duplicates here before continuing
                    W.R =
                        length(sols) > 1 ? unique_points(sols; atol = atol, rtol = rtol) :
                        sols
                end
            end
            update_progress!(progress, W)
        end
    end
end

function regeneration_monodromy_options(
    M::MonodromyOptions,
    W;
    atol = 1e-14,
    rtol = sqrt(eps()),
)

    MonodromyOptions(;
        permutations = false,
        trace_test = true,
        single_loop_per_start_solution = M.single_loop_per_start_solution,
        check_startsolutions = M.check_startsolutions,
        group_actions = M.group_actions,
        loop_finished_callback = M.loop_finished_callback,
        parameter_sampler = M.parameter_sampler,
        equivalence_classes = M.equivalence_classes,
        duplicate_check = M.duplicate_check,
        certification_max_precision = M.certification_max_precision,
        certification_refine_solution = M.certification_refine_solution,
        trace_test_tol = M.trace_test_tol,
        target_solutions_count = Int(floor(1.5 * degree(W))), # in case we miss solutions
        timeout = M.timeout,
        min_solutions = M.min_solutions,
        max_loops_no_progress = M.max_loops_no_progress,
        reuse_loops = M.reuse_loops,
        distance = M.distance,
        # Forward the practical tolerances so monodromy's internal add! uses them
        # rather than the tight per-solution uniqueness_rtol(res) based on ω.
        unique_points_atol = something(M.unique_points_atol, atol),
        unique_points_rtol = something(M.unique_points_rtol, rtol),
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
    warning::Bool = true,
    kwargs...,
)
    progress = cache.progress

    F, f, h = cache.Fᵢ, cache.fᵢ, cache.h # eqs for W and H
    u = cache.u

    P = points(W)
    vars = variables(F)

    # Step 1:
    # we check which points of W are also contained in H
    # the rest is removed from P = points(W) and added to P_next
    # for further processing
    m = .!(is_contained(W, H, cache; threading = threading, kwargs...))
    P_next = manage_initial_points!(P, m)
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
    u_data, d =
        set_up_u_homotopy(W, f, X, h, vars, u; projective = cache.projective, ℓ = cache.ℓ)
    trackers = initialize_u_homotopy_trackers(u_data, cache)


    # the start solutions are the Cartesian product between P_next and the d-th roots of unity.
    roots = [exp(2 * pi * im * k / d) for k = 0:(d-1)]

    # here comes the loop for tracking
    update_progress!(
        progress;
        is_membership_test = false,
        is_solving = true,
        is_monodromy = false,
    )
    if threading
        threaded_intersection!(X, P_next, roots, trackers, u_data, cache, progress, warning)
    else
        serial_intersection!(X, P_next, roots, trackers, u_data, cache, progress, warning)
    end

    nothing
end

function manage_initial_points!(P, m)
    P_next = P[m]
    deleteat!(P, m)
    return P_next
end


function is_contained(X::WitnessPoints, Y::WitnessSet, F, cache; kwargs...)
    progress = cache.progress
    update_progress!(
        progress;
        is_membership_test = true,
        is_solving = false,
        is_monodromy = false,
    )
    # main idea: for every x∈X we take a linear space L with codim(L)=dim(Y) through p and move the points in Y to L. Then, we check if the computed points contain x. If yes, return true, else return false.
    LX = linear_subspace(X)
    LY = linear_subspace(Y)

    if isempty(points(Y)) || isempty(points(X))
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


# The following piece of the code implements the core part of u-regeneration
# It consists of three subhomotopies: Hom1, Hom2, Hom3.
# Hom1 moves the linear equation u = c to a random linear space L1 involving (x,u), where x are the variables of F.
# Hom2 moves u^d - 1 (or u^d - ℓ^d in the projective case) to the new polynomial fᵢ.
# Hom3 moves the linear space L1 to the linear space given by the initial witness sets.
# We first run Hom1. Then we check if the output of Hom1 works as input for Hom2. If this fails, resample L1. Only then track Hom2 and Hom3.

const HOM2_START_CHECK_TARGET = 0.99

is_accepted(T::EndgameTracker, code) = is_success(code) && !T.state.singular

function scaling(p, ℓ_coeffs)
    isnothing(ℓ_coeffs) ? one(ComplexF64) :
    sum(ℓ_coeffs[j] * p[j] for j = 1:length(ℓ_coeffs))
end

function start_solution(P, roots, idx, ℓ_coeffs = nothing)
    l_P = length(P)
    p_idx = rem(idx - 1, l_P) + 1
    r_idx = div(idx - 1, l_P) + 1
    p = copy(P[p_idx])
    p[end] = scaling(p, ℓ_coeffs) * roots[r_idx]
    p
end

function track_intersection!(
    X,
    P,
    roots,
    trackers,
    u_data,
    cache,
    progress,
    threading,
    warning,
)
    l_start = length(P) * length(roots)
    P1 = Vector{Vector{ComplexF64}}(undef, l_start)
    # Vector{Bool} rather than BitVector: the threaded Hom1 sweep writes entries from several threads, and BitVector writes are read-modify-writes of shared chunks.
    accepted = fill(false, l_start)

    ntrials = max(cache.max_trials_u_homotopy, 1)
    for trial = 1:ntrials
        trial > 1 && reset_u_homotopy_trackers!(trackers, u_data, cache)

        fill!(accepted, false)
        # Hom1 and the cheap Hom2 start check are paired: as soon as an
        # intermediate slice produces a start that Hom2 cannot consume, the
        # slice is considered bad and we resample it. On the final trial there
        # is nothing left to resample, so we finish the sweep and keep the
        # starts that did pass the check.
        fail_fast = trial < ntrials
        start_progress_paths!(progress, "Compute start solutions", l_start)
        track_hom1_and_check_hom2_starts!(
            accepted,
            P1,
            P,
            roots,
            trackers[1],
            trackers[2],
            progress,
            threading,
            fail_fast,
            cache.ℓ_coeffs,
        )
        all(accepted) && break
    end

    if warning && any(.!accepted)
        @warn "Some homotopy paths in the u-regeneration intersection step failed. The returned witness set may be incomplete."
    end
    # `accepted` guards against undefined or rejected P1 entries. This matters
    # both after failed trials and when the final trial only validates a subset
    # of the starts.
    start_progress_paths!(progress, "Track paths", l_start)
    track_hom2_hom3!(X, P1, accepted, trackers[2], trackers[3], progress, threading)

    nothing
end

function serial_intersection!(X, P, roots, trackers, u_data, cache, progress, warning)
    track_intersection!(X, P, roots, trackers, u_data, cache, progress, false, warning)
end

function threaded_intersection!(X, P, roots, trackers, u_data, cache, progress, warning)
    track_intersection!(X, P, roots, trackers, u_data, cache, progress, true, warning)
end

function u_homotopy_tracker(H, cache)
    EndgameTracker(
        H;
        tracker_options = cache.tracker_options,
        options = cache.endgame_options,
    )
end

function initialize_u_homotopy_trackers(u_data, cache)
    F, G, L, L2 = u_data
    projective = cache.projective
    L1 = rand_subspace(ambient_dim(L); dim = dim(L), affine = !projective)

    Hom1 = linear_subspace_homotopy(
        F,
        L,
        L1;
        homogeneous = projective,
        intrinsic = cache.intrinsic,
    )
    Hom2 = StraightLineHomotopy(slice(F, L1), slice(G, L1); gamma = cis(2 * pi * rand()))
    projective && (Hom2 = on_affine_chart(Hom2))
    Hom3 = linear_subspace_homotopy(
        G,
        L1,
        L2;
        homogeneous = projective,
        intrinsic = cache.intrinsic,
    )

    [
        u_homotopy_tracker(Hom1, cache),
        u_homotopy_tracker(Hom2, cache),
        u_homotopy_tracker(Hom3, cache),
    ]
end

function reset_u_homotopy_trackers!(trackers, u_data, cache)
    F, G, L, L2 = u_data
    projective = cache.projective
    L1 = rand_subspace(ambient_dim(L); dim = dim(L), affine = !projective)

    target_parameters!(trackers[1], L1)
    Hom2 = StraightLineHomotopy(slice(F, L1), slice(G, L1); gamma = cis(2 * pi * rand()))
    projective && (Hom2 = on_affine_chart(Hom2))
    trackers[2] = u_homotopy_tracker(Hom2, cache)
    start_parameters!(trackers[3], L1)

    trackers
end

function set_up_u_homotopy(W, f, X, h, vars, u; projective::Bool = false, ℓ = 1)

    P, Q = get_num_den(h)
    d = degree(P)
    if projective
        h0 = u^d - ℓ^d
    else
        h0 = u^d - 1
    end
    if degree(Q) > 0
        h0 = h0 / Q
    end
    # we start with the linear space L which does not use pose conditions on u, so that u^d=1
    # we end with the linear space L2 with u=c.
    L = linear_subspace_u(W)
    L2 = linear_subspace(X)

    F = fixed(System([f; h0], variables = vars); compile = false)
    G = fixed(System([f; h], variables = vars); compile = false)

    return (F, G, L, L2), d
end


function track_hom1_and_check_hom2_starts!(
    accepted,
    P1,
    P,
    roots,
    tracker1,
    tracker2,
    progress,
    threading,
    fail_fast,
    ℓ_coeffs = nothing,
)
    if threading
        threaded_track_hom1_and_check_hom2_starts!(
            accepted,
            P1,
            P,
            roots,
            tracker1,
            tracker2,
            progress,
            fail_fast,
            ℓ_coeffs,
        )
    else
        serial_track_hom1_and_check_hom2_starts!(
            accepted,
            P1,
            P,
            roots,
            tracker1,
            tracker2,
            progress,
            fail_fast,
            ℓ_coeffs,
        )
    end
end

function track_hom2_hom3!(X, P1, accepted, tracker2, tracker3, progress, threading)
    if threading
        threaded_track_hom2_hom3!(X, P1, accepted, tracker2, tracker3, progress)
    else
        serial_track_hom2_hom3!(X, P1, accepted, tracker2, tracker3, progress)
    end
end

function serial_track_hom1_and_check_hom2_starts!(
    accepted,
    P1,
    P,
    roots,
    tracker1,
    tracker2,
    progress,
    fail_fast,
    ℓ_coeffs = nothing,
)
    l_start = length(P) * length(roots)

    for i = 1:l_start
        update_progress_paths!(progress, i, l_start)

        p = start_solution(P, roots, i, ℓ_coeffs)
        code1 = track!(tracker1, p, 1)
        if !is_accepted(tracker1, code1)
            accepted[i] = false
            fail_fast && break
            continue
        end

        p1 = solution(tracker1)
        code2 = track!(tracker2.tracker, p1, 1, HOM2_START_CHECK_TARGET)
        if !is_success(code2)
            accepted[i] = false
            fail_fast && break
            continue
        end

        accepted[i] = true
        P1[i] = p1
    end

    nothing
end

function serial_track_hom2_hom3!(X, P1, accepted, tracker2, tracker3, progress)
    l_start = length(P1)

    for i = 1:l_start
        update_progress_paths!(progress, i, l_start)
        accepted[i] || continue

        code2 = track!(tracker2, P1[i], 1)
        if is_accepted(tracker2, code2)
            p2 = solution(tracker2)
            code3 = track!(tracker3, p2, 1)
            if is_accepted(tracker3, code3)
                push!(X, solution(tracker3))
            end
        end
    end

    nothing
end

function threaded_track_hom1_and_check_hom2_starts!(
    accepted,
    P1,
    P,
    roots,
    tracker1,
    tracker2,
    progress,
    fail_fast,
    ℓ_coeffs = nothing,
)

    l_start = length(P) * length(roots)
    nthr = Threads.nthreads()
    trackers1 = [deepcopy(tracker1) for _ = 1:nthr]
    trackers2 = [deepcopy(tracker2.tracker) for _ = 1:nthr]

    progress_lock = ReentrantLock()

    next_idx = Threads.Atomic{Int}(1)
    failed = Threads.Atomic{Bool}(false)
    Threads.@sync begin
        for tid = 1:nthr
            let local_tracker1 = trackers1[tid], local_tracker2 = trackers2[tid]
                Threads.@spawn begin
                    while true
                        fail_fast && failed[] && break

                        idx = Threads.atomic_add!(next_idx, 1)
                        if idx > l_start
                            break
                        end

                        lock(progress_lock) do
                            update_progress_paths!(progress, idx, l_start)
                        end

                        p = start_solution(P, roots, idx, ℓ_coeffs)
                        code1 = track!(local_tracker1, p, 1)
                        if is_accepted(local_tracker1, code1)
                            p1 = solution(local_tracker1)
                            code2 = track!(local_tracker2, p1, 1, HOM2_START_CHECK_TARGET)
                            if is_success(code2)
                                accepted[idx] = true
                                P1[idx] = p1
                            else
                                accepted[idx] = false
                                fail_fast && (failed[] = true)
                            end
                        else
                            accepted[idx] = false
                            fail_fast && (failed[] = true)
                        end
                    end
                end
            end
        end
    end

    nothing
end

function threaded_track_hom2_hom3!(X, P1, accepted, tracker2, tracker3, progress)

    l_start = length(P1)

    nthr = Threads.nthreads()
    trackers2 = [deepcopy(tracker2) for _ = 1:nthr]
    trackers3 = [deepcopy(tracker3) for _ = 1:nthr]

    progress_lock = ReentrantLock()

    next_idx = Threads.Atomic{Int}(1)
    Threads.@sync begin
        for tid = 1:nthr
            let local_tracker2 = trackers2[tid], local_tracker3 = trackers3[tid]
                Threads.@spawn begin
                    while true

                        idx = Threads.atomic_add!(next_idx, 1)
                        if idx > l_start
                            break
                        end

                        lock(progress_lock) do
                            update_progress_paths!(progress, idx, l_start)
                        end
                        accepted[idx] || continue

                        code2 = track!(local_tracker2, P1[idx], 1)
                        if is_accepted(local_tracker2, code2)
                            p2 = solution(local_tracker2)
                            code3 = track!(local_tracker3, p2, 1)
                            if is_accepted(local_tracker3, code3)
                                q = solution(local_tracker3)
                                lock(progress_lock) do
                                    push!(X, q)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    nothing
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
    monodromy_solutions::Int = 0
    monodromy_tracked_loops::Int = 0
    monodromy_generated_loops::Int = 0
    monodromy_no_change::Int = 0
end

function update_progress!(progress::DecomposeProgress; is_monodromy = nothing)
    isnothing(is_monodromy) ? nothing : progress.is_monodromy = is_monodromy
    if is_monodromy === true
        progress.monodromy_solutions = 0
        progress.monodromy_tracked_loops = 0
        progress.monodromy_generated_loops = 0
        progress.monodromy_no_change = 0
    end

    if !isnothing(is_monodromy)
        PM.update!(
            progress.progress_meter,
            progress.step,
            showvalues = showstatus(progress),
        )
    end
end
function update_progress!(
    progress::DecomposeProgress,
    stats::MonodromyStatistics;
    queued::Int,
    solutions::Int,
    finish::Bool = false,
)
    progress.monodromy_solutions = solutions
    progress.monodromy_tracked_loops = stats.tracked_loops[]
    progress.monodromy_generated_loops = stats.generated_loops[]
    progress.monodromy_no_change = loops_no_change(stats, solutions)

    if finish
        PM.update!(
            progress.progress_meter,
            progress.step,
            showvalues = showstatus(progress),
        )
    elseif time() > progress.progress_meter.tlast + progress.progress_meter.dt
        PM.update!(
            progress.progress_meter,
            progress.step,
            showvalues = showstatus(progress),
        )
    end

    nothing
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
            if progress.monodromy_solutions > 0
                push!(
                    text,
                    (
                        "Status",
                        "running monodromy ($(progress.monodromy_solutions) solutions, $(progress.monodromy_generated_loops) loops generated, $(progress.monodromy_no_change) no change)",
                    ),
                )
            else
                push!(text, ("Status", "running monodromy"))
            end
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

## Strategy

The witness points of `W` are grouped into the witness sets of the irreducible components using two kinds of evidence:

* **Monodromy**: moving the points around a loop of slices permutes the points of each irreducible component among themselves.
  Points that are connected by these permutations (an *orbit*) lie on the same component, unless a path jumped between components.
* **The trace test**: the trace of a set of witness points vanishes if and only if the set is a union of complete components (see the comment above [`PointTraces`](@ref)).

The algorithm proceeds as follows.

1. **Completing `W`.** Duplicates are removed and the trace test ([`trace_precheck!`](@ref)) checks whether `W` is complete.
   If not, a monodromy run tries to find the missing points.
2. **Iterations.** As long as there are points that are not certified:
   1. *Traces.* One monodromy loop is run for all remaining points.
      Every new point is tracked once along a fixed pencil of parallel slices `L → L₀₁ → L₁` and its contribution to the trace is cached ([`update_traces!`](@ref)).
      The displacement is the same for all points of `W`, so the trace of every set of points is a sum of cached values ([`cached_trace`](@ref)) and no further path tracking is necessary.
   2. *Orbits.* The permutations of the loop are merged with the orbits of the previous iterations ([`get_orbits_from_monodromy_permutations`](@ref)).
      A link between a small-trace point (whose own trace vanishes) and another point is only used after it was seen in two different loops, so that a single path jump cannot join components.
   3. *Splitting.* Orbits whose trace vanishes but that contain a smaller set of at least 2 points with vanishing trace are split ([`split_nonminimal_orbits`](@ref)).
      Such an orbit contains several components, e.g., after a path jumped.
   4. *Single orbits.* Every orbit that satisfies the certification policy [`is_certifiable`](@ref) is certified.
      If a point of an orbit has no cached trace, its trace is computed by a separate monodromy run instead (fallback).
   5. *Unions of orbits.* A component can consist of several orbits that the loops have not connected yet.
      The combinatorial trace test ([`certify_orbit_unions!`](@ref)) certifies the smallest unions of the remaining orbits with vanishing trace.
   6. The certified points are removed, the orbits are re-indexed, and the iteration stops once the degrees of the certified components add up to the number of points.
3. Points that could not be certified after `max_iters` iterations are returned in witness sets with `is_irreducible = nothing`.

## Certification policy

All certifications (steps 2.4 and 2.5) use [`is_certifiable`](@ref):
* the trace of the set must vanish,
* sets with at least 2 points must be trace-minimal: no subset with at least 2 points (and at least 2 points outside of it) may have a vanishing trace,
* a single point is only certified as a component of degree 1 after `DEGREE1_ITERATIONS` iterations in which all orbits consist of single points and no new points or new links of small-trace points were found (*degree 1 rule*).
  The trace of a single point is only weak evidence, since a single point on a component of degree > 1 can have a small trace.
  It is therefore never used to reject or split larger sets; instead, links of small-trace points need to be confirmed by a second loop (step 2.2).

Phantom points (new points of a monodromy run that coincide with points of an already certified component) are removed from the orbits before certification.
"""
function decompose_with_monodromy!(
    W,
    show_monodromy_progress,
    options,
    max_iters,
    warning,
    progress,
    seed;
    intrinsic::Union{Nothing,Bool} = nothing,
    threading::Bool = Threads.nthreads() > 1,
)
    update_progress_dim!(progress, dim(W))

    P₀ = points(W)
    L = linear_subspace(W)
    G = system(W)
    n = ambient_dim(L)

    # tolerances for comparing points
    atol = something(options.unique_points_atol, 1e-14)
    rtol = something(options.unique_points_rtol, 1e-8)

    # first check that all points in P are unique
    certified_up = UniquePoints(first(P₀), 1; distance = InfNorm())
    P = Vector{eltype(P₀)}()
    for (i, pᵢ) in enumerate(P₀)
        _, new_point = add!(certified_up, pᵢ, i; atol = atol, rtol = rtol)
        if new_point
            push!(P, pᵢ)
        end
    end
    empty!(certified_up) # we reuse this later

    decomposition = Vector{WitnessSet}()

    if dim(L) < n
        update_progress!(progress; is_monodromy = true)

        MS =
            MonodromySolver(G, L; compile = false, options = options, intrinsic = intrinsic)
        # Cache of per-point traces for this witness set.
        # The random displacement for the trace test is chosen here, once, and used for all points of this witness set, so that traces of different points can be added.
        point_traces = PointTraces(L, options)
        # How often a link between a small-trace point and another point was observed in a monodromy loop, see `get_orbits_from_monodromy_permutations`.
        small_trace_links = SmallTraceLinks()
        initial_points = check_start_solutions(MS, P, L)
        # complete the witness set; the permutations of this call are not used
        if trace_precheck!(MS, solution.(initial_points), L, threading)
            non_complete_points = solution.(initial_points)
        else
            res = monodromy_solve(
                MS,
                solution.(initial_points),
                L,
                seed;
                threading = threading,
                show_progress = show_monodromy_progress,
                progress = show_monodromy_progress ? nothing : progress,
            )
            non_complete_points = indexed_solutions(res)
            if warning && (trace(res) > options.trace_test_tol)
                @warn "Trying to decompose non-complete set of witness points for codimension $(dim(L)) (trace test failed). Will try to compute the missing points. The output will contain all components, for which the trace test was successfull."
            end
        end
        update_progress!(progress; is_monodromy = false)
        update_progress_npts!(progress, length(non_complete_points))

        iter = 0
        allow_degree1_iter = 0
        d = length(non_complete_points) # the total degree (i.e., number of points on all irreducible components)
        non_complete_orbits = Vector{Set{Int}}()
        trace_certified_all = false

        while !isempty(non_complete_points)
            iter += 1
            if iter > max_iters
                break
            end


            update_progress!(progress; is_monodromy = true)
            res = monodromy_solve(
                MS,
                non_complete_points,
                L,
                seed;
                threading = threading,
                show_progress = show_monodromy_progress,
                progress = show_monodromy_progress ? nothing : progress,
            )
            update_progress!(progress; is_monodromy = false)

            prev_count = length(non_complete_points)
            updated_points = indexed_solutions(res)
            d += length(updated_points) - length(non_complete_points) # update total degree in case we found new points.
            non_complete_points = updated_points

            # Identify phantom points: new points that jumped into an already-certified component.
            # Only check when new points were found and certified_up is non-empty.
            phantom_indices = Set{Int}()
            if length(updated_points) > prev_count && length(certified_up) > 0
                # reset degree1 check if we found new points 
                allow_degree1_iter = 0
                # now compute phantom points
                for (i, p) in pairs(updated_points)
                    rad = max(atol, norm(p, Inf) * rtol)
                    if !isnothing(search_in_radius(certified_up, p, rad))
                        push!(phantom_indices, i)
                    end
                end
                d -= length(phantom_indices)
            end
            ℓ = length(non_complete_points) # once for a later check
            k = length(non_complete_points) # and once for counting progress
            update_progress_npts!(progress, k)

            # If monodromy deduplicated starting solutions, the cached orbit indices are
            # stale (they reference the old count). Reset to avoid a BoundsError
            # also reset the degree 1 counter in this case
            if length(non_complete_points) < prev_count
                non_complete_orbits = Vector{Set{Int}}()
                allow_degree1_iter = 0
            end

            # Track every point that is new in this iteration once along the fixed pencil L → L₀₁ → L₁.
            # Points from earlier iterations are already cached.
            update_traces!(point_traces, non_complete_points, MS, threading)

            # Orbits from the permutations of this loop and the orbits of the previous iterations.
            # Links involving small-trace points (potential degree 1 points) are only used once they are confirmed by a second loop.
            orbits, new_links = get_orbits_from_monodromy_permutations(
                permutations(res),
                non_complete_orbits,
                non_complete_points,
                small_trace_points(
                    point_traces,
                    non_complete_points,
                    options.trace_test_tol,
                ),
                small_trace_links,
            )
            # Something happened to the potential degree 1 points: restart the waiting period of the degree 1 rule.
            new_links > 0 && (allow_degree1_iter = 0)

            # An orbit can contain points of several components if a monodromy path jumped between components (e.g. close to an intersection point).
            # Split such orbits along subsets with vanishing trace.
            orbits = split_nonminimal_orbits(
                orbits,
                non_complete_points,
                point_traces,
                phantom_indices,
                options.trace_test_tol,
            )

            # Degree 1 rule (see `is_certifiable`): single points are only certified after DEGREE1_ITERATIONS iterations in which all orbits consist of single points, i.e., when there is nothing else to classify.
            # The counter is reset when new points are found or a new link of a small-trace point is observed.
            if all(o -> length(o) == 1, orbits)
                allow_degree1_iter += 1
            end
            complete_orbits = Vector{Set{Int}}()

            for orbit in orbits
                update_progress!(progress; is_monodromy = true)

                # Strip phantom indices from the orbit before running inner monodromy.
                clean_orbit =
                    isempty(phantom_indices) ? orbit : setdiff(orbit, phantom_indices)
                if isempty(clean_orbit)
                    push!(complete_orbits, orbit)
                    k -= length(orbit)
                    update_progress_npts!(progress, k)
                    update_progress!(progress; is_monodromy = false)
                    continue
                end

                P_orbit = non_complete_points[sort!(collect(clean_orbit))]
                # The trace of the orbit is computed from the cached per-point traces, no path tracking is needed.
                # Only if a point has no cached trace (one of its trace paths failed) we fall back to a monodromy run for this orbit, which computes the trace with a new random loop.
                orbit_trace = cached_trace(point_traces, P_orbit)
                if isnothing(orbit_trace)
                    res_orbit = monodromy_solve(
                        MS,
                        P_orbit,
                        L,
                        seed;
                        threading = threading,
                        show_progress = show_monodromy_progress,
                        progress = show_monodromy_progress ? nothing : progress,
                    )
                    orbit_trace = trace(res_orbit)
                    P_orbit = indexed_solutions(res_orbit)
                end
                update_progress!(progress; is_monodromy = false)

                # certify the orbit if it satisfies the certification policy, see `is_certifiable` (this includes the degree 1 rule)
                if is_certifiable(
                    point_traces,
                    P_orbit,
                    orbit_trace,
                    options.trace_test_tol;
                    single_point = length(clean_orbit) == 1,
                    allow_single_point = allow_degree1_iter ≥ DEGREE1_ITERATIONS ||
                                         iter ≥ max_iters - 1,
                )
                    P_certified = P_orbit
                    W_new = WitnessSet(G, L, P_certified; is_irreducible = true)

                    if length(P_certified) > length(clean_orbit)
                        # Inner monodromy found new points beyond the starting orbit.
                        # Match them against non_complete_points to detect orbit merging and path-jumping phantoms.
                        matched = setdiff(
                            Set(
                                matching_indices(
                                    non_complete_points,
                                    P_certified;
                                    atol = atol,
                                    rtol = rtol,
                                ),
                            ),
                            phantom_indices,
                        )
                        # If not all P_certified points are accounted for in non_complete_points, the result is unreliable (phantom or genuine new point) — skip.
                        if length(matched) < length(P_certified)
                            continue
                        end
                        complete_orbit = matched
                    else
                        complete_orbit = clean_orbit
                    end

                    push!(complete_orbits, complete_orbit)
                    k -= length(complete_orbit)
                    update_progress_npts!(progress, k)

                    push!(decomposition, W_new)
                    for p in solutions(W_new)
                        add!(certified_up, p, 1; atol = atol, rtol = rtol)
                    end
                    d += length(P_certified) - length(complete_orbit)
                    update_progress!(progress, W_new)
                end
            end

            # Combinatorial trace test: certify unions of the remaining orbits whose trace vanishes.
            # This finds a component as soon as all its points are in some orbits, even if the loops have not merged them into a single orbit yet.
            # See `certify_orbit_unions!`.
            certify_orbit_unions!(
                complete_orbits,
                decomposition,
                certified_up,
                orbits,
                phantom_indices,
                non_complete_points,
                point_traces,
                G,
                L,
                options.trace_test_tol;
                atol = atol,
                rtol = rtol,
            ) do orbit, W_new
                k -= length(orbit)
                update_progress_npts!(progress, k)
                update_progress!(progress, W_new)
            end

            # Check if we are done
            if sum(degree, decomposition; init = 0) == d
                trace_certified_all = true
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
                append!(complete_orbit_indices, orbit)
            end
            append!(complete_orbit_indices, phantom_indices)
            sort!(complete_orbit_indices)
            unique!(complete_orbit_indices)

            deleteat!(non_complete_points, complete_orbit_indices)


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

            # Skip removed indices (complete orbits and phantoms) when shifting.
            shift_orbit(orbit) = Set(
                orbit_indices_mapping[i] for i in orbit if haskey(orbit_indices_mapping, i)
            )
            # 1. shift existing non complete orbits to new mapping
            non_complete_orbits =
                shift_orbit.(
                    merge_sets(
                        [
                            # Remove all orbits that overlap with a complete orbit or phantom indices
                            filter(non_complete_orbits) do o
                                all(co -> isdisjoint(co, o), complete_orbits) &&
                                    isdisjoint(o, phantom_indices)
                            end
                            # A complete orbit can differ from the orbit it came from (phantoms removed, matched points, unions of orbits), so remove all orbits that overlap a complete orbit.
                            filter(orbits) do o
                                all(co -> isdisjoint(co, o), complete_orbits)
                            end
                        ],
                    ),
                )
        end

        if !trace_certified_all && !isempty(non_complete_points)
            covered_indices = Set{Int}()
            for orbit in non_complete_orbits
                union!(covered_indices, orbit)
            end

            orbits_to_return = copy(non_complete_orbits)
            for i in eachindex(non_complete_points)
                if !(i in covered_indices)
                    push!(orbits_to_return, Set([i]))
                end
            end

            for orbit in orbits_to_return
                isempty(orbit) && continue
                P_orbit = non_complete_points[sort!(collect(orbit))]
                W_new = WitnessSet(G, L, P_orbit; is_irreducible = nothing)

                push!(decomposition, W_new)
                update_progress!(progress, W_new)
            end

            if warning
                @warn "Some witness sets could not be certified by the trace test and are returned with is_irreducible = nothing. They will not be displayed in a numerical irreducible decomposition."
            end
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

##
## Trace test with cached per-point traces
##
## Background.
## Let W be a witness set on the slice L and let L(s) = L + s⋅v, s ∈ ℂ, be a family of parallel slices.
## Every witness point x moves along a path x(s) with x(0) = x.
## For a set S of witness points let s_S(s) = Σ_{x ∈ S} x(s).
## Then:

## S is a union of complete irreducible components  ⟺  s_S(s) is affine-linear in s
##
## (for generic L and v).
## We evaluate s_S at s = 0, 1, 2, i.e. on the slices L, L₀₁ = L + v and L₁ = L + 2v, and measure how far the three points s_S(0), s_S(1), s_S(2) are from lying on a line.
## This is the trace test of `monodromy_solve`, see `trace_colinearity`.
##
## Since s_S is a sum over the points of S, each point only has to be tracked once along L → L₀₁ → L₁.
## Afterwards the trace of every subset S is a sum of cached values and no path tracking is needed.
## This is used in `decompose_with_monodromy!` to test orbits and unions of orbits (`certify_orbit_unions!`).

"""
    trace_precheck!(MS::MonodromySolver, X, L::LinearSubspace, threading::Bool)

Check with the trace test whether the witness points `X` on the slice `L` are complete, without running a monodromy loop.
For this the points are only tracked along the part of a loop that collects the trace, `L → L₀₁ → L₁`.
Returns `true` if the trace test succeeds.
Then no monodromy is necessary to complete the witness points.
The state of `MS` is reset afterwards, so `monodromy_solve(MS, ...)` can be called directly.
Returns `false` if some path failed, the trace test fails, or no point of `X` is a valid start solution.
"""
function trace_precheck!(
    MS::MonodromySolver,
    X::AbstractVector{<:AbstractVector{<:Number}},
    L::LinearSubspace,
    threading::Bool,
)
    reset!() = begin
        MS.statistics = MonodromyStatistics()
        empty!(MS.unique_points)
        duplicate_check_is_certified(MS) && empty!(MS.certified_solutions)
        reset_trace!(MS)
        reset_loops!(MS)
    end
    reset!()
    cL = convert(LinearSubspace{ComplexF64}, L)
    results =
        MS.options.check_startsolutions ? check_start_solutions(MS, X, cL) :
        initialize_start_solutions(MS, X)
    isempty(results) && return false
    add_loop!(MS, cL)
    ℓ = loop(MS, 1)
    ok = Threads.Atomic{Bool}(true)
    nthr = threading ? min(Threads.nthreads(), length(results)) : 1
    if nthr > 1
        for i = length(MS.trackers)+1:nthr
            push!(MS.trackers, deepcopy(MS.trackers[1]))
        end
        next_idx = Threads.Atomic{Int}(1)
        Threads.@sync for tid = 1:nthr
            let tracker = MS.trackers[tid].tracker
                Threads.@spawn while ok[]
                    i = Threads.atomic_add!(next_idx, 1)
                    i > length(results) && break
                    x₁ = track_trace_legs!(tracker, results[i], ℓ, MS.trace, MS.trace_lock)
                    isnothing(x₁) && (ok[] = false)
                end
            end
        end
    else
        tracker = MS.trackers[1].tracker
        for r in results
            x₁ = track_trace_legs!(tracker, r, ℓ, MS.trace, MS.trace_lock)
            isnothing(x₁) && (ok[] = false; break)
        end
    end
    success = ok[] && trace_colinearity(MS) < MS.options.trace_test_tol
    reset!()
    success
end

"""
    PointTraces(L::LinearSubspace, options::MonodromyOptions)

Cache for the trace test in [`decompose`](@ref) for one witness set with slice `L`.

On construction a random displacement `v` is drawn once (via `MonodromyLoop(L, options)`, which sets `L₀₁ = L + v` and `L₁ = L + 2v`).
It is never changed afterwards: the field `loop` is immutable and every trace in the cache is computed with this same `loop`.
This is essential, since traces of different points may only be added if they are computed for the same displacement.

* `data[x]` is an `n × 3` matrix with the positions of the witness point `x` on `L`, `L₀₁` and `L₁` (columns 1, 2, 3), or `nothing` if one of the two paths failed.

A `PointTraces` must only be used for points on the slice `L` of the witness set it was created for.
In `decompose_with_monodromy!` it is a local variable, created once per witness set, so different witness sets never share a cache.
"""
struct PointTraces
    loop::MonodromyLoop{LinearSubspace{ComplexF64}}
    data::Dict{Vector{ComplexF64},Union{Nothing,Matrix{ComplexF64}}}
end
PointTraces(L::LinearSubspace, options::MonodromyOptions) = PointTraces(
    MonodromyLoop(L, options),
    Dict{Vector{ComplexF64},Union{Nothing,Matrix{ComplexF64}}}(),
)

"""
    point_trace(tracker::Tracker, x, loop::MonodromyLoop)

Track `x` from `loop.p` over `loop.p₀₁` to `loop.p₁` and return the `n × 3` matrix of the three positions, or `nothing` if a path failed.
"""
function point_trace(tracker::Tracker, x::Vector{ComplexF64}, loop::MonodromyLoop)
    T = zeros(ComplexF64, length(x), 3)
    x₁ = track_trace_legs!(tracker, x, loop, T, nothing)
    isnothing(x₁) ? nothing : T
end

"""
    update_traces!(PT::PointTraces, P, MS::MonodromySolver, threading::Bool)

Compute and cache the traces of all points in `P` that are not cached yet.
The paths are tracked with the trackers of `MS` (one per thread), always along `PT.loop`.
"""
function update_traces!(PT::PointTraces, P, MS::MonodromySolver, threading::Bool)
    # copy the points, since they are used as keys of the cache
    new_points = [copy(p) for p in P if !haskey(PT.data, p)]
    isempty(new_points) && return PT
    traces = Vector{Union{Nothing,Matrix{ComplexF64}}}(undef, length(new_points))
    nthr = threading ? min(Threads.nthreads(), length(new_points)) : 1
    if nthr > 1
        for i = length(MS.trackers)+1:nthr
            push!(MS.trackers, deepcopy(MS.trackers[1]))
        end
        next_idx = Threads.Atomic{Int}(1)
        Threads.@sync for tid = 1:nthr
            let tracker = MS.trackers[tid].tracker
                Threads.@spawn while true
                    i = Threads.atomic_add!(next_idx, 1)
                    i > length(new_points) && break
                    traces[i] = point_trace(tracker, new_points[i], PT.loop)
                end
            end
        end
    else
        tracker = MS.trackers[1].tracker
        for (i, p) in enumerate(new_points)
            traces[i] = point_trace(tracker, p, PT.loop)
        end
    end
    for (p, T) in zip(new_points, traces)
        PT.data[p] = T
    end
    PT
end

"""
    cached_trace(PT::PointTraces, S)

The trace of the set of witness points `S`, computed from the cached traces in `PT`: with `s₀, s₀₁, s₁` the sums of the positions of the points of `S` on `L, L₀₁, L₁`, this is `σ₃ / σ₁` for the singular values of the `(n+1) × 3` matrix `[s₀ s₀₁ s₁; 1 1 1]`.
It is (numerically) zero if and only if the three sums lie on a line, i.e., if `S` is a union of complete components.
This is exactly the value that `monodromy_solve` reports as `trace`.

Returns `nothing` if `S` is empty or a point of `S` has no cached trace (not computed yet, or a path failed).
"""
function cached_trace(PT::PointTraces, S)
    isempty(S) && return nothing
    n = length(first(S))
    M = zeros(ComplexF64, n + 1, 3)
    M[end, :] .= 1
    for p in S
        T = get(PT.data, p, nothing)
        isnothing(T) && return nothing
        for j = 1:3, i = 1:n
            M[i, j] += T[i, j]
        end
    end
    σ = LA.svdvals(M)
    σ[3] / σ[1]
end

##
## Certification policy
##
## `decompose_with_monodromy!` certifies a set of witness points as one irreducible component in three places: single orbits, unions of orbits (`certify_orbit_unions!`), and the points of a fallback monodromy run.
## All of them use `is_certifiable`, so that the rules below are applied consistently.
##
## The evidence comes from monodromy (an orbit lies on one component, unless a path jumped) and from the trace test.
## The trace of a single point is only *weak* evidence: a single point on a component of degree > 1 can have a small trace (e.g. where the component is almost linear).
## The trace of a set with at least 2 points is *strong* evidence.
## Accordingly, a single point's trace is never used to reject or split a larger set.
##
## Points with a small trace of their own (*small-trace points*, potential degree 1 points) are ambiguous: they are a component of degree 1, or a point on a larger component.
## The degree 1 rule resolves this with monodromy:
## * a link (a path of a loop) between a small-trace point and another point is only used to form orbits after it was observed in 2 different loops (`get_orbits_from_monodromy_permutations`).
##   A single path jump therefore does not join two components, while the points of a larger component are connected repeatedly;
## * a small-trace point that is not linked to other points is certified as a component of degree 1 only after `DEGREE1_ITERATIONS` iterations without new points and without new links of small-trace points, see `is_certifiable`.

# Number of iterations (in which all orbits consist of single points) before a single point is certified as a component of degree 1.
const DEGREE1_ITERATIONS = 15

# Sets with at most this many points are checked for trace-minimality using all subsets, larger sets only using subsets with 2 points.
const MAX_MINIMALITY_CHECK_POINTS = 12

"""
    is_certifiable(PT::PointTraces, S, trace_S, tol;
        single_point = length(S) == 1, allow_single_point = false)

The certification policy of [`decompose`](@ref): decide whether the witness points `S` with trace `trace_S` are certified as the witness points of one irreducible component.

1. The trace must vanish: `trace_S < tol` (`trace_S = nothing` means that the trace is not known, then `S` is not certified).
2. Degree 1 rule: if the evidence is a single point (`single_point`), it is only certified if `allow_single_point`.
   In `decompose_with_monodromy!` this is the case after `DEGREE1_ITERATIONS` iterations in which all orbits are single points and no new points and no new links of small-trace points were found, i.e., when the monodromy loops had enough chances to connect the point to other points.
   Links of small-trace points need a confirmation by a second loop before a small-trace point joins a larger orbit (see [`get_orbits_from_monodromy_permutations`](@ref)), so a small-trace point cannot escape this rule by a single path jump.
3. Trace-minimality: if `S` has at least 2 points, no proper subset of `S` with at least 2 points may have a vanishing trace (see [`passing_proper_subset`](@ref)).
   A union of complete components has a vanishing trace, so a passing proper subset shows that `S` consists of several components.
   Subsets with a single point are not used here, since their trace is only weak evidence (see 2.).

Rule 3 cannot separate components of degree 1 from other components (for example two lines joined by a path jump), since this would require the weak evidence of single points.
These cases are handled by rule 2 together with the confirmation of links in [`get_orbits_from_monodromy_permutations`](@ref).
"""
function is_certifiable(
    PT::PointTraces,
    S,
    trace_S,
    tol;
    single_point::Bool = length(S) == 1,
    allow_single_point::Bool = false,
)
    (isnothing(trace_S) || trace_S ≥ tol) && return false
    single_point && return allow_single_point
    isnothing(passing_proper_subset(PT, S, tol))
end

"""
    passing_proper_subset(PT::PointTraces, S, tol)

Return the positions (in `S`) of a proper subset of the witness points `S` with at least 2 points (and at least 2 points outside of it) whose cached trace is smaller than `tol`, or `nothing` if there is none.
Subsets are tested by increasing size, so a smallest one is returned.

This is only meaningful if the trace of `S` vanishes.
Then the trace of the complement of a passing subset vanishes as well, so only subsets with at most `length(S) ÷ 2` points are tested.
Points without a cached trace never belong to a passing subset.
If `S` has more than `MAX_MINIMALITY_CHECK_POINTS` points, only subsets with 2 points are tested.
"""
function passing_proper_subset(PT::PointTraces, S, tol)
    k = length(S)
    k ≤ 3 && return nothing
    if k ≤ MAX_MINIMALITY_CHECK_POINTS
        masks = filter(m -> 2 ≤ count_ones(m) ≤ k ÷ 2, collect(1:(1<<k)-2))
        sort!(masks; by = count_ones)
        subsets = ([j for j = 1:k if (m >> (j - 1)) & 1 == 1] for m in masks)
    else
        subsets = ([i, j] for i = 1:k for j = i+1:k)
    end
    for sub in subsets
        t = cached_trace(PT, S[sub])
        !isnothing(t) && t < tol && return sub
    end
    nothing
end

"""
    split_nonminimal_orbits(orbits, P, PT::PointTraces, phantom_indices, tol)

Split every orbit whose trace vanishes but that is not trace-minimal (it has a passing proper subset, see [`passing_proper_subset`](@ref)) into that subset and the rest, and repeat this for the parts.
Such an orbit contains points of several components, e.g., after a path jumped between components.
The parts can then be certified separately.
Orbits with a phantom point, without cached traces or with a non-vanishing trace are kept as they are.
"""
function split_nonminimal_orbits(orbits, P, PT::PointTraces, phantom_indices, tol)
    out = Vector{Set{Int}}()
    queue = collect(Set{Int}, orbits)
    while !isempty(queue)
        orbit = pop!(queue)
        if length(orbit) > 3 && isdisjoint(orbit, phantom_indices)
            idx = sort!(collect(orbit))
            t = cached_trace(PT, P[idx])
            if !isnothing(t) && t < tol
                sub = passing_proper_subset(PT, P[idx], tol)
                if !isnothing(sub)
                    A = Set(idx[sub])
                    push!(queue, A, setdiff(orbit, A))
                    continue
                end
            end
        end
        push!(out, orbit)
    end
    sort!(out; by = length, rev = true)
end

"""
    small_trace_points(PT::PointTraces, P, tol)

The indices of the points of `P` whose own cached trace is smaller than `tol`.
These are the *potential degree 1 points*: the point is either a component of degree 1 by itself or a point on a larger component whose trace happens to be small.
Points without a cached trace are not included.
"""
function small_trace_points(PT::PointTraces, P, tol)
    small = Set{Int}()
    for (i, p) in pairs(P)
        t = cached_trace(PT, (p,))
        !isnothing(t) && t < tol && push!(small, i)
    end
    small
end

"""
    SmallTraceLinks()

For each pair of witness points `{p, q}` where at least one of them is a small-trace point, the number of monodromy loops in which a path connected `p` and `q`.
The pairs are stored by the coordinates of the points, so the counts are independent of the re-indexing of the points between iterations.
"""
const SmallTraceLinks = Dict{Set{Vector{ComplexF64}},Int}

# Number of different loops that have to connect a small-trace point to another point before this link is used to form an orbit.
const SMALL_TRACE_LINK_CONFIRMATIONS = 2

"""
    get_orbits_from_monodromy_permutations(perms, initial_orbits, P, small, links::SmallTraceLinks)

Compute the orbits of the witness points `P` from the permutations `perms` (one column per loop, `0` marks a failed path) and the orbits of the previous iterations `initial_orbits`.
Returns the orbits (sorted by decreasing size) and the number of links involving small-trace points that were observed for the first time.

Two points are in the same orbit if they are connected by links, where a link is either
* a remembered orbit from `initial_orbits`, or
* a path of a loop from `i` to `perms[i]`.

A link of a path that involves a small-trace point (`i ∈ small` or `perms[i] ∈ small`) is only used if it was observed in at least `SMALL_TRACE_LINK_CONFIRMATIONS` different loops (counted in `links`, which is updated).
A single observation can be a path jump, and a small-trace point in a larger orbit escapes the degree 1 rule (see [`is_certifiable`](@ref)).
Requiring a confirmation by an independent loop separates two components that were joined by a single jump (e.g. two lines), while the points of a component of higher degree are connected repeatedly.

Since the returned orbits are remembered as `initial_orbits` in the next iteration, a small-trace point is only remembered in an orbit with other points via confirmed links.
"""
function get_orbits_from_monodromy_permutations(
    perms::AbstractMatrix,
    initial_orbits,
    P,
    small::Set{Int},
    links::SmallTraceLinks,
)
    parent = Dict{Int,Int}()
    rank = Dict{Int,Int}()
    add!(i) = haskey(parent, i) || (parent[i] = i; rank[i] = 0)

    for orbit in initial_orbits
        isempty(orbit) && continue
        r = first(orbit)
        add!(r)
        for i in orbit
            add!(i)
            merge_sets_union(parent, rank, i, r)
        end
    end

    new_links = 0
    for col in eachcol(perms)
        # every pair is counted at most once per loop
        observed = Set{Set{Vector{ComplexF64}}}()
        for (i, j) in enumerate(col)
            j == 0 && continue
            add!(i)
            add!(j)
            i == j && continue
            if i in small || j in small
                key = Set([P[i], P[j]])
                if !(key in observed)
                    push!(observed, key)
                    links[key] = get(links, key, 0) + 1
                    links[key] == 1 && (new_links += 1)
                end
                links[key] ≥ SMALL_TRACE_LINK_CONFIRMATIONS || continue
            end
            merge_sets_union(parent, rank, i, j)
        end
    end

    groups = Dict{Int,Set{Int}}()
    for i in keys(parent)
        push!(get!(groups, merge_sets_find(parent, i), Set{Int}()), i)
    end
    sort!(collect(values(groups)); by = length, rev = true), new_links
end

# The combinatorial trace test tests up to 2^MAX_COMBINATORIAL_ORBITS unions of orbits.
const MAX_COMBINATORIAL_ORBITS = 10

"""
    certify_orbit_unions!(on_certified, complete_orbits, decomposition, certified_up,
        orbits, phantom_indices, P, point_traces, G, L, tol; atol, rtol)

Combinatorial trace test: find unions of monodromy orbits that are exactly the witness points of one irreducible component, and certify them.

### Why this works
1. A monodromy loop moves a witness point only to points on the *same* irreducible component, so every orbit lies in a single component.
   But a component can consist of several orbits, as long as the loops have not connected them yet.
2. The trace of a set of witness points vanishes if and only if the set is a union of complete components (see the comment above `PointTraces`).

Let `U` be a union of orbits whose trace vanishes, such that no union of fewer of these orbits has a vanishing trace.
If `U` contained points of two components `A` and `B`, then `A ∩ U` would be a union of fewer orbits (by 1.) with vanishing trace (by 2.), a contradiction.
So `U` is the set of witness points of one irreducible component.

### Algorithm
1. Candidates are the orbits (without the indices in `phantom_indices`) that
   * do not overlap an orbit in `complete_orbits` (already certified),
   * have a cached trace for every point (see [`cached_trace`](@ref)),
   * do **not** have a vanishing trace themselves.
     The last condition is needed for the minimality argument above: an orbit with vanishing trace is a candidate for a component on its own (for example a single point on a line) and is handled by the single-orbit certification, including the degree 1 rule (see [`is_certifiable`](@ref)).
     If two such orbits were candidates here, their union would have a vanishing trace as well and could be certified as one component.
2. If there are fewer than 2 or more than `MAX_COMBINATORIAL_ORBITS` candidates nothing is done.
3. All unions of at least 2 candidates are enumerated as bit masks, sorted by the number of orbits they contain (all pairs, then all triples, ...).
   A union that shares an orbit with a union certified before is skipped.
   If the union satisfies the certification policy [`is_certifiable`](@ref) (vanishing trace and trace-minimality), it is certified: it is added to `complete_orbits`, a new irreducible [`WitnessSet`](@ref) with system `G` and slice `L` is added to `decomposition`, its points are added to `certified_up`, and `on_certified(indices, W_new)` is called.

Testing the unions by increasing size guarantees the minimality required above.
No path tracking happens here; all traces come from `point_traces`.

### Arguments
* `P`: the witness points; orbits are sets of indices into `P`.
* `atol`, `rtol`: tolerances for adding points to `certified_up`.
"""
function certify_orbit_unions!(
    on_certified,
    complete_orbits::Vector{Set{Int}},
    decomposition,
    certified_up,
    orbits,
    phantom_indices,
    P,
    point_traces::PointTraces,
    G,
    L,
    tol;
    atol,
    rtol,
)
    covered = Set{Int}()
    for o in complete_orbits
        union!(covered, o)
    end

    # 1. candidate orbits
    remaining = Vector{Vector{Int}}()
    for orbit in orbits
        o = sort!(collect(setdiff(orbit, phantom_indices)))
        (isempty(o) || !isdisjoint(o, covered)) && continue
        t = cached_trace(point_traces, P[o])
        # no cached trace, or a component on its own
        (isnothing(t) || t < tol) && continue
        push!(remaining, o)
    end

    # 2. size limit
    m = length(remaining)
    (2 ≤ m ≤ MAX_COMBINATORIAL_ORBITS) || return nothing

    # 3. bit j-1 of `mask` set ⟺ the union contains the orbit remaining[j]
    masks =
        sort!(filter(mask -> count_ones(mask) ≥ 2, collect(1:(1<<m)-1)); by = count_ones)
    used = 0 # orbits that are part of a certified union
    for mask in masks
        mask & used == 0 || continue
        idx = Int[]
        for j = 1:m
            (mask >> (j - 1)) & 1 == 1 && append!(idx, remaining[j])
        end
        sort!(idx)
        t = cached_trace(point_traces, P[idx])
        # The union is minimal among unions of candidate orbits (by the order of the masks); `is_certifiable` additionally checks it against all subsets of its points with at least 2 points.
        if is_certifiable(point_traces, P[idx], t, tol)
            used |= mask
            W_new = WitnessSet(G, L, P[idx]; is_irreducible = true)
            push!(complete_orbits, Set(idx))
            push!(decomposition, W_new)
            for p in P[idx]
                add!(certified_up, p, 1; atol = atol, rtol = rtol)
            end
            on_certified(idx, W_new)
        end
    end
    nothing
end

# Helper function to merge two sets
function merge_sets_find(parent, i)
    root = i
    while parent[root] != root
        root = parent[root]
    end
    while parent[i] != root
        next = parent[i]
        parent[i] = root
        i = next
    end
    root
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

# check if we have found points that are already classifiedxw
# find indices of newly detected points
function matching_indices(points, certified_points; atol, rtol)
    indices = Int[]
    matched = falses(length(certified_points))
    for (i, p) in pairs(points)
        rad = max(atol, norm(p, Inf) * rtol)
        for (j, q) in pairs(certified_points)
            matched[j] && continue
            if distance(p, q, InfNorm()) ≤ rad
                push!(indices, i)
                matched[j] = true
                break
            end
        end
    end
    indices
end


function decompose_with_monodromy_options(
    M::MonodromyOptions;
    atol = 1e-14,
    rtol = sqrt(eps()),
)

    # we need two copies of the options.
    # one with trace test
    # one without to run monodromy populating the permutation matrix
    options = MonodromyOptions(;
        permutations = true,
        trace_test = true,
        single_loop_per_start_solution = true,
        check_startsolutions = false,
        group_actions = M.group_actions,
        loop_finished_callback = M.loop_finished_callback,
        parameter_sampler = M.parameter_sampler,
        equivalence_classes = M.equivalence_classes,
        duplicate_check = M.duplicate_check,
        certification_max_precision = M.certification_max_precision,
        certification_refine_solution = M.certification_refine_solution,
        trace_test_tol = M.trace_test_tol,
        target_solutions_count = M.target_solutions_count,
        timeout = M.timeout,
        min_solutions = M.min_solutions,
        max_loops_no_progress = M.max_loops_no_progress,
        reuse_loops = M.reuse_loops,
        distance = M.distance,
        unique_points_atol = something(M.unique_points_atol, atol),
        unique_points_rtol = something(M.unique_points_rtol, rtol),
    )

    options
end


"""
    decompose(Ws::Vector{WitnessPoints}; options...) 

This function decomposes a [`WitnessSet`](@ref) or a vector of [`WitnessSet`](@ref) into irreducible components.

### Options
* `show_progress = true`: indicate whether a progress bar should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) should be displayed. If `false`, minimal info about the monodromy computations are still displayed in the progress bar of `decompose`.
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_iters = 5000`: maximal number of iterations for the decomposition step.
* `warning = true`: if `true`, prints a warning when the [`trace_test`](@ref) fails.
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
    Ws::Vector{WP};
    show_progress::Bool = true,
    show_monodromy_progress::Bool = false,
    monodromy_options::MonodromyOptions = MonodromyOptions(),
    max_iters::Int = 5000,
    intrinsic::Union{Nothing,Bool} = nothing,
    warning::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    seed = nothing,
    atol = 1e-14,
    rtol = sqrt(eps()),
) where {WP<:WitnessSet}

    if isnothing(seed)
        seed = rand(UInt32)
    end
    Random.seed!(seed)

    sort!(Ws; by = dim, rev = true)
    options = decompose_with_monodromy_options(monodromy_options; atol = atol, rtol = rtol)
    out = Vector{WitnessSet}()

    if isempty(Ws)
        return out
    end

    c = length(Ws)
    n = ambient_dim(linear_subspace(Ws[1]))

    if show_progress
        progress_meter = PM.ProgressUnknown(
            dt = 0.1,
            desc = "Decomposing $c witness sets",
            enabled = true,
            spinner = true,
        )
        progress_meter.tlast += 0.3
        progress = DecomposeProgress(progress_meter = progress_meter, n_witness_sets = c)
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
                intrinsic = intrinsic,
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
struct NumericalIrreducibleDecomposition{T<:WitnessSet}
    witness_sets::Dict{Int,Vector{T}}
    seed::Union{Nothing,UInt32}
end
NumericalIrreducibleDecomposition(Ws::Vector{T}) where {T<:WitnessSet} =
    NumericalIrreducibleDecomposition(Ws, nothing)
function NumericalIrreducibleDecomposition(Ws::Vector{T}, seed) where {T<:WitnessSet}
    D = Dict{Int,Vector{T}}()
    for W in Ws
        k = dim(W)
        push!(get!(D, k, T[]), W)
    end
    NumericalIrreducibleDecomposition(D, seed)
end

"""
    witness_sets(N::NumericalIrreducibleDecomposition;
                dims::Union{Vector{Int},Nothing} = nothing,
                only_irreducible::Bool = true)

Returns the witness sets in `N`. 
`dims` specifies the dimensions that should be considered.
"""
function witness_sets(
    N::NumericalIrreducibleDecomposition{T};
    dims::Union{Vector{Int},Nothing} = nothing,
    dim::Union{Int,Nothing} = nothing,
    only_irreducible::Bool = true,
) where {T}
    W = N.witness_sets
    out = Dict{Int,Vector{T}}()
    selected_dims = isnothing(dim) ? dims : [dim]
    selected_dims = isnothing(selected_dims) ? keys(W) : selected_dims

    for k in selected_dims
        if haskey(W, k)
            Ws = only_irreducible ? filter(_is_irreducible, W[k]) : W[k]
            if !isempty(Ws)
                out[k] = Ws
            end
        end
    end

    out
end
witness_sets(
    N::NumericalIrreducibleDecomposition{T},
    dim::Int;
    only_irreducible::Bool = true,
) where {T} = witness_sets(N; dim = dim, only_irreducible = only_irreducible)
seed(N::NumericalIrreducibleDecomposition{T}) where {T} = N.seed

"""
    ncomponents(N::NumericalIrreducibleDecomposition;
                dims::Union{Vector{Int},Nothing} = nothing,
                only_irreducible::Bool = true)

Returns the total number of components in `N`. 
`dims` specifies the dimensions that should be considered.
"""
function ncomponents(
    N::NumericalIrreducibleDecomposition{T};
    dims::Union{Vector{Int},Nothing} = nothing,
    dim::Union{Int,Nothing} = nothing,
    only_irreducible::Bool = true,
) where {T}
    dims = isnothing(dim) ? dims : [dim]
    D = witness_sets(N; dims = dims, only_irreducible = only_irreducible)
    if isempty(D)
        return 0
    else
        return sum(length(last(Ws)) for Ws in D)
    end
end
ncomponents(
    N::NumericalIrreducibleDecomposition{T},
    dim::Int;
    only_irreducible::Bool = true,
) where {T} = ncomponents(N; dim = dim, only_irreducible = only_irreducible)
n_components(N; dims = nothing, dim = nothing, only_irreducible::Bool = true) =
    ncomponents(N; dims = dims, dim = dim, only_irreducible = only_irreducible)
n_components(N, dim; only_irreducible::Bool = true) =
    ncomponents(N; dim = dim, only_irreducible = only_irreducible)

"""

    degrees(N::NumericalIrreducibleDecomposition;
            dims::Union{Vector{Int},Nothing} = nothing,
            only_irreducible::Bool = true)

Returns the degrees of the components in `N`.
`dims` specifies the dimensions that should be considered.

"""
function ModelKit.degrees(
    N::NumericalIrreducibleDecomposition{T};
    dims::Union{Vector{Int},Nothing} = nothing,
    dim::Union{Int,Nothing} = nothing,
    only_irreducible::Bool = true,
) where {T}
    dims = isnothing(dim) ? dims : [dim]
    D = witness_sets(N; dims = dims, only_irreducible = only_irreducible)
    out = Dict{Int,Vector{Int}}()
    for key in keys(D)
        out[key] = degree.(D[key])
    end

    out
end
ModelKit.degrees(
    N::NumericalIrreducibleDecomposition{T},
    dim::Int;
    only_irreducible::Bool = true,
) where {T} = degrees(N; dim = dim, only_irreducible = only_irreducible)

function max_dim(
    N::NumericalIrreducibleDecomposition{T};
    only_irreducible::Bool = true,
) where {T}
    D = witness_sets(N; only_irreducible = only_irreducible)
    k = keys(D)
    if !isempty(k)
        maximum(k)
    else
        -1
    end
end

function Base.show(io::IO, N::NumericalIrreducibleDecomposition{T}) where {T}
    D = witness_sets(N)
    total = ncomponents(N)
    s = total == 1 ? "component" : "components"
    header = "Numerical irreducible decomposition with $total $s"
    println(io, header)
    println(io, "="^(length(header)))
    mdim = max_dim(N)
    if mdim >= 0
        for d = mdim:-1:0
            if haskey(D, d)
                ℓ = count(_is_irreducible, D[d])
                if ℓ > 0
                    println(io, "• $ℓ component(s) of dimension $d.")
                end
            end
        end

        println(io, "\n degree table of components:")
        degree_table(io, N)
    end
end
function degree_table(io, N::NumericalIrreducibleDecomposition{T}) where {T}
    D = degrees(N)
    k = collect(keys(D))
    sort!(k, rev = true)
    n = length(k)

    headers = ["dimension", "degrees of components"]
    data = Matrix{Union{Int,String}}(undef, n, 2)

    for (i, key) in enumerate(k)
        data[i, 1] = key
        components = D[key]
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
* `sorted`: if `true` (default), the polynomials in `F`` will be sorted by degree in increasing order (when `F` is a system of rational functions, we first compute witness sets for each entry of `F` and then sort `F` according to their degrees). If `false, the polynomials in `F` will not be sorted. The third option is `:randomized`, which multiplies `F` by a random matrix.
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `monodromy_options_for_regeneration`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref) in [`regeneration`](@ref).
* `monodromy_options_for_decompose`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref) in [`decompose`](@ref).
* `show_monodromy_progress = false`: if `true`, sets `show_monodromy_for_regeneration_progress` and `show_monodromy_for_decompose_progress` to `true`. If `false`, minimal info about the monodromy computations are still displayed in the progress bar of each substep.
* `show_monodromy_for_regeneration_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) in [`regeneration`](@ref) should be displayed. 
* `show_monodromy_for_decompose_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) in [`decompose`](@ref) should be displayed.
* `max_iters = 5000`: maximal number of iterations for the decomposition step.
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).` This option is used for [`regeneration`](@ref).
* `warning = true`: if `true`, prints warnings when the [`trace_test`](@ref) or an intermediate u-regeneration intersection step fails.
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
    show_progress::Bool = true,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_accuracy = 1e-10,
    ),
    monodromy_options_for_regeneration = MonodromyOptions(),
    monodromy_options_for_decompose::MonodromyOptions = MonodromyOptions(),
    show_monodromy_progress::Bool = false,
    show_monodromy_for_regeneration_progress::Bool = false,
    show_monodromy_for_decompose_progress::Bool = false,
    max_iters::Int = 5000,
    sorted::Union{Bool,Symbol} = true,
    max_codim::Union{Int,Nothing} = nothing,
    max_trials_u_homotopy::Int = 5,
    intrinsic_for_regeneration::Union{Nothing,Bool} = true,
    intrinsic_for_decompose::Union{Nothing,Bool} = nothing,
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
        show_progress = show_progress,
        sorted = sorted,
        max_codim = max_codim,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
        monodromy_options = monodromy_options_for_regeneration,
        max_trials_u_homotopy = max_trials_u_homotopy,
        intrinsic = intrinsic_for_regeneration,
        show_monodromy_progress = show_monodromy_for_regeneration_progress,
        threading = threading,
        seed = nothing,
        atol = atol,
        rtol = rtol,
        warning = warning,
        kwargs...,
    )
    if isnothing(Ws)
        return nothing
    end
    dec = decompose(
        Ws;
        show_progress = show_progress,
        monodromy_options = monodromy_options_for_decompose,
        max_iters = max_iters,
        intrinsic = intrinsic_for_decompose,
        show_monodromy_progress = show_monodromy_for_decompose_progress,
        threading = threading,
        warning = warning,
        seed = seed,
        atol = atol,
        rtol = rtol,
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
    path_label::String = "Track paths"
    current_task::Int = 0
    ntasks::Int = 0
    current_path::Int = 0
    npaths::Int = 0
    monodromy_solutions::Int = 0
    monodromy_tracked_loops::Int = 0
    monodromy_generated_loops::Int = 0
    monodromy_no_change::Int = 0
end
IntersectProgress(progress_meter::PM.ProgressUnknown) =
    IntersectProgress(progress_meter = progress_meter)
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
    if is_monodromy === true
        progress.monodromy_solutions = 0
        progress.monodromy_tracked_loops = 0
        progress.monodromy_generated_loops = 0
        progress.monodromy_no_change = 0
    end
    PM.update!(progress.progress_meter, showvalues = showvalues(progress))
end
function update_progress!(
    progress::IntersectProgress,
    stats::MonodromyStatistics;
    queued::Int,
    solutions::Int,
    finish::Bool = false,
)
    progress.monodromy_solutions = solutions
    progress.monodromy_tracked_loops = stats.tracked_loops[]
    progress.monodromy_generated_loops = stats.generated_loops[]
    progress.monodromy_no_change = loops_no_change(stats, solutions)

    if finish
        PM.update!(progress.progress_meter, showvalues = showvalues(progress))
    elseif time() > progress.progress_meter.tlast + progress.progress_meter.dt
        PM.update!(progress.progress_meter, showvalues = showvalues(progress))
    end

    nothing
end
update_progress!(progress::IntersectProgress, W) = nothing
function update_progress_monodromy_witness!(progress::IntersectProgress, W)
    progress.monodromy_solutions = degree(W)
    progress.monodromy_tracked_loops = 0
    progress.monodromy_generated_loops = 0
    progress.monodromy_no_change = 0
    PM.update!(progress.progress_meter, showvalues = showvalues(progress))
end
function update_progress_tasks!(progress::IntersectProgress, i::Int, m::Int)
    progress.current_task = i
    progress.ntasks = m
    PM.update!(progress.progress_meter, showvalues = showvalues(progress))
end
function update_progress_paths!(progress::IntersectProgress, i::Int, m::Int)
    progress.current_path = i
    progress.npaths = m
    PM.update!(progress.progress_meter, showvalues = showvalues(progress))
end
function start_progress_paths!(progress::IntersectProgress, label::String, m::Int)
    progress.path_label = label
    progress.current_path = 0
    progress.npaths = m
    PM.update!(progress.progress_meter, showvalues = showvalues(progress))
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
        if progress.is_solving && progress.npaths > 0
            push!(
                text,
                (progress.path_label, "$(progress.current_path) / $(progress.npaths)"),
            )
        elseif progress.is_monodromy
            push!(
                text,
                (
                    "Fill up points",
                    "$(progress.monodromy_solutions) solutions, $(progress.monodromy_tracked_loops) loops tracked ($(progress.monodromy_generated_loops) generated, $(progress.monodromy_no_change) no change)",
                ),
            )
        elseif progress.is_membership_test
            push!(
                text,
                ("Membership test", "$(progress.current_task) / $(progress.ntasks)"),
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

    fᵢ::Vector{Expression}
    Fᵢ::Sys
    h::Expression
    u::Variable
    projective::Bool
    ℓ::Union{Nothing,Expression}
    ℓ_coeffs::Union{Nothing,Vector{ComplexF64}}

    max_trials_u_homotopy::Int
    intrinsic::Union{Nothing,Bool}

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{IntersectProgress,Nothing}
end

function IntersectCache(
    u,
    f,
    F,
    h,
    projective,
    ℓ,
    ℓ_coeffs,
    max_trials_u_homotopy,
    intrinsic,
    EO,
    TO,
    progress,
)
    m, N = size(F)
    A = zeros(ComplexF64, N - 1, N)
    b = zeros(ComplexF64, N - 1)
    x0 = zeros(ComplexF64, N)
    y0 = zeros(ComplexF64, m)
    y = zeros(ComplexF64, m)

    IntersectCache(
        A,
        b,
        x0,
        y0,
        y,
        f,
        F,
        h,
        u,
        projective,
        ℓ,
        ℓ_coeffs,
        max_trials_u_homotopy,
        intrinsic,
        EO,
        TO,
        progress,
    )
end


"""
    intersect(W::WitnessSet, H::WitnessSet)

This intersects the witness sets `W` and `H`, where `H` is a hypersurface defined by single polynomial `f`. 

### Options

* `show_progress = true`: indicate whether a progress bar should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress bar of [`monodromy_solve`](@ref) should be displayed. If `false`, minimal info about the monodromy computations are still displayed in the progress bar of `intersect`.
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_trials_u_homotopy = 5`: maximal number of random subspaces tried until an intermediate u-regeneration intersection step succeeds.
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).`
* `warning = true`: if `true`, prints a warning when an intermediate u-regeneration intersection step fails after all trials.
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)

### Example

The following example computes witness sets for a union of two circles.

```julia-repl
julia> @var x y z
julia> H = witness_set(x^2 + y^2 + z^2 - 1)
julia> W = witness_set([x^2 + y^2 - 4; (x-1)^2 + (z-1)^2 - 1], codim = 2)
julia> intersect(W, H)
Witness set for dimension 0 of degree 8
```    

"""
function Base.intersect(W::WitnessSet, H::WitnessSet; kwargs...)
    out = _intersect(W, H; kwargs...)
    if length(out) == 1
        return first(out)
    else
        return out
    end
end

"""
    intersect(W::WitnessSet, f::Expression)

First computes a witness set `H` for `f` and then runs  `intersect(W, H)`.
"""
function Base.intersect(W::WitnessSet, f::Expression; kwargs...)

    vars = variables(system(W))
    @assert all(v -> v ∈ vars, variables(f)) "Witness sets must be in the same ambient space."

    F = System([f], variables = vars)
    H = witness_set(F)
    intersect(W, H; kwargs...)
end

function _intersect(
    W::WitnessSet,
    H::WitnessSet;
    show_progress::Bool = true,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_cond = 1e12,
    ),
    monodromy_options::MonodromyOptions = MonodromyOptions(),
    max_trials_u_homotopy::Int = 5,
    intrinsic::Union{Nothing,Bool} = true, # always intrinsic linear subspace homotopies
    show_monodromy_progress::Bool = false,
    threading = Threads.nthreads() > 1,
    atol = 1e-14,
    rtol = sqrt(eps()),
    warning::Bool = true,
    kwargs...,
)
    @assert size(system(H), 1) == 1 "The second argument must be defined by a single polynomial."
    @assert variables(system(W)) == variables(system(H)) "Witness sets must be in the same ambient space."
    W.projective == H.projective ||
        throw(ArgumentError("Witness sets must both be affine or projective."))

    # progress bar
    if show_progress
        progress_meter = PM.ProgressUnknown(
            dt = 0.4,
            desc = "Intersecting...",
            enabled = true,
            spinner = true,
        )
        progress_meter.tlast += 0.3
        progress = IntersectProgress(progress_meter)
    else
        progress = nothing
    end


    # transform W and H so that they use the additional variable u
    n = ambient_dim(W.L)
    @unique_var u
    vars = variables(system(W))
    vars_u = [vars; u]
    projective = W.projective
    W₁, W₂, Hᵤ, f, F, h = prepare_for_u_homotopy(H, W, vars, vars_u, projective)


    # cache
    ℓ_coeffs = projective ? randn(ComplexF64, n) : nothing
    ℓ = projective ? sum(ℓ_coeffs .* vars) : nothing
    cache = IntersectCache(
        u,
        f,
        F,
        h,
        projective,
        ℓ,
        ℓ_coeffs,
        max_trials_u_homotopy,
        intrinsic,
        endgame_options,
        tracker_options,
        progress,
    )

    # intersect
    intersect_with_hypersurface!(
        W₁,
        Hᵤ,
        W₂,
        cache;
        threading = threading,
        atol = atol,
        rtol = rtol,
        warning = warning,
        kwargs...,
    )
    update_Fᵢ!(cache, [f; h], vars_u)
    fill_up!(
        [W₁; W₂],
        monodromy_options,
        cache,
        show_monodromy_progress,
        threading;
        atol = atol,
        rtol = rtol,
    )

    # return data 
    G = fixed(System([f; h], variables = vars); compile = false)
    P1, L1 = u_transform(W₁)
    out = [WitnessSet(G, L1, P1; projective = projective)]
    if !isnothing(W₂)
        P2, L2 = u_transform(W₂)
        push!(out, WitnessSet(G, L2, P2; projective = projective))

    end
    filter!(X -> degree(X) > 0, out)

    finish_progress!(progress)
    out
end



function prepare_for_u_homotopy(H, W, vars, vars_u, projective)
    # prepare polynomials 
    FH0 = deepcopy(System(system(H)))
    FW0 = deepcopy(System(system(W)))
    h = subs(expressions(FH0), variables(FH0) => vars)
    f = subs(expressions(FW0), variables(FW0) => vars)
    FH = System(h, variables = vars_u)
    FW = System(f, variables = vars_u)

    # prepare flags
    LW = linear_subspace(W)
    LH = linear_subspace(H)
    flagW = get_flag(1:2, LW) # returns L, Lᵤ
    flagH = get_flag(1:1, LH)
    cW, cH = get_c(flagW), get_c(flagH)

    # prepare WitnessPoints
    # in u-coordiantes the points in W have to be transformed as x -> [x; c]
    W₁ = WitnessPoints(flagW[1][1], flagW[1][2], map(x -> [x; cW], solutions(W)))
    if dim(W) == 0
        W₂ = nothing
    else
        W₂ = WitnessPoints(flagW[2][1], flagW[2][2], Vector{Vector{ComplexF64}}())
    end
    Hᵤ = WitnessSet(
        fixed(FH; compile = false),
        flagH[1][1],
        map(x -> [x; cH], solutions(H));
        projective = projective,
    )

    W₁, W₂, Hᵤ, f, fixed(FW; compile = false), first(h)
end
get_c(flag) = extrinsic((flag[1][1])).b[1] # the first entry of b is the right-hand side of "u=c"
