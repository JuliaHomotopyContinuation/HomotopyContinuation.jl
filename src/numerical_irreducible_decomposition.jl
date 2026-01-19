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
mutable struct WitnessPoints{
    Sub1<:AbstractSubspace,
    Sub2<:AbstractSubspace,
    R<:Vector{ComplexF64},
}
    L::Sub1
    Lᵤ::Sub2
    R::Vector{R}
end
linear_subspace(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = W.L
linear_subspace_u(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = W.Lᵤ
points(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = W.R
points(W::WitnessSet{A,B,C}) where {A,B,C} = W.R
codim(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = dim(linear_subspace(W))
dim(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = codim(linear_subspace(W))
ModelKit.degree(W::WitnessPoints{A,B,Vector{ComplexF64}}) where {A,B} = length(points(W))
push!(W::WitnessPoints{A,B,Vector{ComplexF64}}, R::Vector{T}) where {A,B,T<:Number} =
    push!(W.R, R)

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
    is_removing_points::Bool = false
    is_computing_hypersurfaces::Bool = true
    current_task::Int = 0
    ntasks::Int = 0
    current_hypersurface::Int = 1
    nhypersurfaces::Int
    current_dim::Int = 0
    next_dim::Int = 0
end
WitnessSetsProgress(n::Int, nhypersurfaces::Int, progress_meter::PM.ProgressUnknown) =
    WitnessSetsProgress(
        ambient_dim = n,
        nhypersurfaces = nhypersurfaces,
        progress_meter = progress_meter,
    )
update_progress!(progress::Nothing; 
                is_computing_hypersurfaces::Bool = true,
                is_solving::Bool = false,
                is_removing_points::Bool = false) = nothing
function update_progress!(progress::WitnessSetsProgress; 
                is_computing_hypersurfaces::Bool = false,
                is_solving::Bool = false,
                is_removing_points::Bool = false) 
    progress.is_computing_hypersurfaces = is_computing_hypersurfaces
    progress.is_solving = is_solving
    progress.is_removing_points = is_removing_points
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
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
update_progress_tasks!(progress::Nothing, i::Int, m::Int) = nothing
function update_progress_tasks!(progress::WitnessSetsProgress, i::Int, m::Int)
    progress.current_task = i
    progress.ntasks = m
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress_dim!(progress::Nothing, d::Int) = nothing
function update_progress_dim!(progress::WitnessSetsProgress, d::Int)
    progress.current_dim = d
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress!(progress::Nothing, W::Union{WitnessPoints,WitnessSet}) = nothing
function update_progress!(progress::WitnessSetsProgress, W::Union{WitnessPoints,WitnessSet})
    progress.degrees[codim(W)] = length(points(W))
    PM.update!(
        progress.progress_meter,
        progress.current_hypersurface,
        showvalues = showvalues(progress),
    )
end
update_progress!(progress::Union{Nothing,WitnessSetsProgress}, W::Nothing) = nothing
finish_progress!(progress::Nothing) = nothing
function finish_progress!(progress::WitnessSetsProgress)
    progress.is_solving = false
    progress.is_removing_points = false
    PM.finish!(progress.progress_meter, showvalues = showvalues(progress))
end
function showvalues(progress::WitnessSetsProgress)
    if progress.is_computing_hypersurfaces
        text = [("Computing hypersurfaces", "")]
    else
        text = [(
            "Intersect with hypersurface",
            "$(progress.current_hypersurface) / $(progress.nhypersurfaces)",
        )]
        if progress.is_removing_points
            push!(
                text,
                ("Remove junk points from dim. $(progress.current_dim)", "$(progress.current_task)/$(progress.ntasks)"),
            )
        elseif progress.is_solving
            push!(
                text,
                ("Tracking paths", "$(progress.current_task)/$(progress.ntasks)"),
            )
        end
                   
        push!(text, ("Current number of witness points", ""))
        for c = 1:progress.current_hypersurface
            d = progress.ambient_dim - c
            deg = get(progress.degrees, c, nothing)
            if !isnothing(deg) && deg > 0
                push!(text, ("Dimension $d", "$(deg)"))
            end
        end
    end

    text
end

"""
    RegenerationCache

A cache for [`regeneration`](@ref).
"""
mutable struct RegenerationCache{Sys<:AbstractSystem}
    As::Vector
    bs::Vector
    x0::Vector
    y0::Vector
    y::Vector
    U::UniquePoints

    Fᵢ::Sys
    u::Variable
    n::Int
    i::Int
    codim::Int

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{WitnessSetsProgress,Nothing}
end

function RegenerationCache(Fᵢ, u, n, codim, EO, TO, progress)
    As = [zeros(ComplexF64, n + 1 - i, n + 1) for i = 0:codim]
    bs = [zeros(ComplexF64, n + 1 - i) for i = 0:codim]
    x0 = zeros(ComplexF64, n + 1)
    y0 = zeros(ComplexF64, length(Fᵢ))
    y = zeros(ComplexF64, length(Fᵢ))
    U = UniquePoints(x0, 0)

    RegenerationCache(As, bs, x0, y0, y, U, Fᵢ, u, n, 0, codim, EO, TO, progress)
end
function update_Fᵢ!(cache, Fᵢ)
    cache.Fᵢ = Fᵢ
    cache.y0 = zeros(ComplexF64, length(Fᵢ))
    cache.y = zeros(ComplexF64, length(Fᵢ))

    nothing
end
update_i!(cache, i) = cache.i = i
function update_x0!(x0)
    for i = 1:length(x0)
        x0[i] = randn(ComplexF64)
    end
    LA.normalize!(x0)
    x0
end

"""
    regeneration(F::System; options...) 

This solves ``F=0`` equation-by-equations and returns a [`WitnessSet`](@ref) for every dimension without decomposing them into irreducible components (witness sets that are not decomposed are also called witness supersets).

The implementation is based on the algorithm [u-regeneration](https://arxiv.org/abs/2206.02869) by Duff, Leykin and Rodriguez. 

### Options

* `sorted = true`: the polynomials in F will be sorted by degree in decreasing order. 
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `show_progress = true`: indicate whether the progress of the computation should be displayed.
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
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
    threading::Bool = true,
    seed = nothing,
)

    !isnothing(seed) && Random.seed!(seed)

    # the algorithm is u-regeneration as proposed 
    # by Duff, Leykin and Rodriguez in https://arxiv.org/abs/2206.02869

    vars = variables(F)
    if sorted
        f = sort(expressions(F), by = ModelKit.degree, rev = true)
    else
        f = expressions(F)
    end

    n = length(vars) # ambient dimension
    c = length(f) # we can have witness sets of codimesion at most min(c,n)
    expected_max_codim = min(c, n)
    if !isnothing(max_codim) && max_codim < expected_max_codim
        # if max_codim is smaller than the expected codimension we must compute witness points for one more codimension, so that we can remove spurious points
        codim = max_codim + 1
    else
        codim = expected_max_codim
    end

    if show_progress
        progress = WitnessSetsProgress(
            n,
            c,
            PM.ProgressUnknown(
                dt = 1.0,
                desc = "Computing witness sets...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end
    

    # u-regeneration adds another variable u to F
    @unique_var u
    push!(vars, u)

    # initialize the linear equations for witness sets
    A, b, Aᵤ, bᵤ = initialize_linear_equations(n)

    # prepare c witness sets for the output
    # internally we represent a witness superset by WitnessPoints
    # the i-th witness superset out[i] is for codimension i
    out = initialize_witness_sets(codim, n, A, b, Aᵤ, bᵤ)


    # we compute witness (super)sets for the hypersurfaces f[1]=0,...,f[c]=0.
    # it is covenient to use the WitnessSet wrapper here, because this also keeps track of the equation
    # as a linear subspace we take the linear subspace for out[1], that sets u=0.
    update_progress!(progress; is_computing_hypersurfaces = true)
    H = initialize_hypersurfaces(f, vars, linear_subspace(out[1]); threading = threading)
    if any(H .== nothing)
        return nothing
    end

    # Initialize a cache
    Fᵢ = fixed(System(f[1:1], variables = vars), compile = false)
    cache = RegenerationCache(Fᵢ, u, n, codim, endgame_options, tracker_options, progress)

    # now comes the core loop of the algorithm.
    # we start with the first hypersurface f[1]=0 and take its witness superset H[1]
    # then, we for i in 2:c we iteratively intersect all current witness sets with f[i]=0
    update_progress!(progress;
        is_computing_hypersurfaces = false,
        is_solving = true)
    begin
        for i = 1:c
            update_progress_hypersurface!(progress, i)
            if i == 1
                # the first step: take the witness superset H[1] as initial witness superset
                begin
                    P = solutions(H[1])
                    X = out[1]
                    for p in P
                        push!(X, p)
                    end
                    update_progress!(progress, X)
                end
            else
                begin
                    update_progress!(progress;
                                    is_solving = true,
                                    is_removing_points = false)
                    update_i!(cache, i)

                    # for all W in out we intersect W with H[i]
                    intersect_all!(out, H, cache, threading)

                    # update Fᵢ
                    Fᵢ = fixed(System(f[1:i], variables = vars), compile = false)
                    update_Fᵢ!(cache, Fᵢ)

                    update_progress!(progress;
                                    is_solving = false,
                                    is_removing_points = true)
                    # after the first loop that takes care of intersecting with Hᵢ
                    # we now check if we have added points that are already contained in 
                    # witness sets of higher dimension.
                    # we only need to do this for witness sets of codimensions 0<k<n.
                    remove_points_all!(out, cache; threading = threading)

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
    if !isempty(out)
        return map(out) do W
            P, L = u_transform(W)
            WitnessSet(fixed(F, compile = false), L, P)
        end
    else
        return Vector{WitnessSet}()
    end

end


function initialize_linear_equations(n)

    A₀ = randn(ComplexF64, n - 1, n)
    b₀ = randn(ComplexF64, n - 1)
    svd = LA.svd(A₀) # we orthogonalize A
    Aᵥ = svd.Vt

    # u regenerates operates on 2 types of linear equations Ax=b
    # type 1 does not use u, so we set the last column to zero
    Aᵤ = [Aᵥ zeros(n - 1)]
    bᵤ = inv.(svd.S) .* (svd.U' * b₀)
    # type 2 sets the linear equation u=0. We add this as the first equation
    A = [zeros(1, n) 1.0; Aᵥ zeros(n - 1)]
    b = [0; bᵤ]

    (A, b, Aᵤ, bᵤ)
end
function initialize_witness_sets(codim, n, A, b, Aᵤ, bᵤ)
    out = Vector{WitnessPoints}(undef, codim)
    for i = 1:codim
        j = i + 1
        # type 1 linear equation for out[i] takes the last i rows of Aᵤ
        Eᵤ = ExtrinsicDescription(Aᵤ[i:end, :], bᵤ[i:end]; orthonormal = true)
        Lᵤ = LinearSubspace(Eᵤ)
        # type 2 linear equation for out[i] takes the last i-1 rows and the first row of A
        E = ExtrinsicDescription(A[[1; j:n], :], b[[1; j:n]]; orthonormal = true)
        L = LinearSubspace(E)

        out[i] = WitnessPoints(L, Lᵤ, Vector{Vector{ComplexF64}}())
    end

    out
end
function initialize_hypersurfaces(f::Vector{Expression}, vars, L; threading::Bool = true)
    c = length(f)
    out = Vector{WitnessSet}(undef, c)

    for i = 1:c
        h = fixed(System([f[i]], variables = vars), compile = false)
        res = solve(
            h,
            target_subspace = L;
            start_system = :total_degree,
            show_progress = false,
            threading = threading,
        )
        if isnothing(res)
            return nothing
        end
        S = solutions(res, only_nonsingular = true)

        out[i] = WitnessSet(h, L, S)
    end

    out
end

function intersect_all!(out, H, cache, threading)

    i = cache.i
    codim = cache.codim
    progress = cache.progress

    # the i-th hypersurface
    Hᵢ = H[i]

    # we enumerate reversely, so that we can add points to witness sets that we have already
    # taken care of; i.e., if we intersect Wₖ∩Hᵢ below in the intersect_with_hypersurface function,
    # we have already intersected Wₖ₊₁ ∩ Hᵢ.
    update_progress!(progress;
                    is_solving = true,
                    is_removing_points = false)
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
            intersect_with_hypersurface!(Wₖ, Hᵢ, Wₖ₊₁, cache; threading = threading)
            update_progress!(progress, Wₖ)
            update_progress!(progress, Wₖ₊₁)
        end
    end
end


function remove_points_all!(out, cache; threading::Bool = true)

    i = cache.i
    progress = cache.progress
    current_task = 0

    E = enumerate(out)
    for (k, W) in E # k = codim(W) for W in Ws
        update_progress_dim!(progress, dim(W))
        if k > 1 && k <= i && degree(W) > 0
            for j = 1:(k-1)
                current_task +=1
                update_progress_tasks!(progress, current_task, k-1)

                X = out[j]
                if degree(X) > 0
                    remove_points!(W, X, cache; threading = threading)
                    update_progress!(progress, W)
                end
            end
        end
        current_task = 0
        update_progress!(progress, W)
    end
end


"""
    intersect_with_hypersurface!

This is the core routine of the regeneration algorithm. It intersects a set of [`WitnessPoints`](@ref) with a hypersurface.
"""
function intersect_with_hypersurface!(W, H, X, cache; threading::Bool = true)

    F = cache.Fᵢ
    progress = cache.progress
    u = cache.u

    P = points(W)
    f = (System(F).expressions) # equations for W
    G = system(H) # H is the hypersurface
    g = System(G).expressions # equations for H

    # Step 1:
    # we check which points of W are also contained in H
    # the rest is removed from P = points(W) and added to P_next
    # for further processing
    m = .!(is_contained(W, H, cache; threading = threading))
    P_next = manage_initial_points!(P, m, W, progress)

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
    update_progress!(progress, W)
    return P_next
end

function serial_intersection!(X, start, tracker, progress)

    l_start = length(start)
    for (i, s) in enumerate(start)
        p = s[1]
        p[end] = s[2] # the last entry of s[1] is zero. we replace it with a d-th root of unity.

        res = track(tracker, p, 1)
        if is_success(res) && is_finite(res) && is_nonsingular(res)
            new = copy(tracker.state.solution)
            push!(X, new)
            update_progress!(progress, X)
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
            let local_tracker = trackers[tid]
                Threads.@spawn begin
                    while true
                        idx = Threads.atomic_add!(next_idx, 1)
                        if idx > length(start_vec)
                            break
                        end

                        s = start_vec[idx]
                        p = copy(s[1])
                        p[end] = s[2] # the last entry of s[1] is zero. we replace it with a d-th root of unity.

                        res = track(local_tracker, p, 1)

                        if is_success(res) && is_finite(res) && is_nonsingular(res)
                            new = copy(local_tracker.state.solution)
                            lock(progress_lock) do
                                push!(X, new)
                                update_progress!(progress, X)
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

    d = ModelKit.degree(H)
    γ = exp(2 * pi * im * rand()) # gamma trick
    g0 = γ * (u^d - 1)

    # we start with the linear space L which does not use pose conditions on u, so that u^d=1
    # we end with the linear space K with u=0.
    L = linear_subspace_u(W)
    K = linear_subspace(X)

    F₀ = slice(System([f; g0], variables = vars), L; compile = false)
    G₀ = slice(System([f; g], variables = vars), K; compile = false)
    Hom = StraightLineHomotopy(F₀, G₀)

    return Hom, d
end

"""
    is_contained(X, Y, F, cache)

Returns a boolean vector indicating whether the points of X are contained in (Y, F).
"""
function is_contained(X, Y, F, cache; threading::Bool = true)

    tracker_options = cache.tracker_options
    endgame_options = cache.endgame_options


    # main idea: for every x∈X we take a linear space L with codim(L)=dim(Y) through p and move the points in Y to L. Then, we check if the computed points contain x. If yes, return true, else return false.
    LX = linear_subspace(X)
    LY = linear_subspace(Y)
    k = codim(LY) - codim(LX) # k≥0 iff dim X ≤ dim Y

    if k < 0 || length(points(Y)) == 0 || length(points(X)) == 0
        # if dim X > dim Y return only false
        out = falses(length(points(X)))
    else
        # setup 
        U = cache.U
        empty!(U)

        Hom = linear_subspace_homotopy(F, LY, LY; intrinsic = true)
        tracker = EndgameTracker(
            Hom;
            tracker_options = tracker_options,
            options = endgame_options,
        )
        set_up_linear_spaces!(cache, LX, LY)

        if threading
            out = threaded_x_in_Y(X, Y, F, tracker, cache)
        else
            out = serial_x_in_Y(X, Y, F, tracker, cache)
        end
    end

    out
end
function is_contained(
    V::WitnessPoints,
    W::WitnessSet,
    cache::RegenerationCache;
    threading::Bool = true,
)
    is_contained(V, W, system(W), cache; threading = threading)
end

function remove_points!(
    W::WitnessPoints,
    V::WitnessPoints,
    cache::RegenerationCache;
    threading::Bool = true,
)
    m = is_contained(W, V, cache.Fᵢ, cache; threading = threading)
    deleteat!(W.R, m)

    nothing
end

function set_up_linear_spaces!(cache, LX, LY)

    n = ambient_dim(LY)
    cX = codim(LX)
    dY = dim(LY)
    k = codim(LY) - codim(LX)

    # to compute linear spaces through the points in X we first set up
    # a matrix-vector pair of the correct size cY×n
    A = cache.As[dY+1]
    b = cache.bs[dY+1]

    # since we have used a subset of the equations for Y also for X,
    # we can reuse them
    AX = extrinsic(LX).A
    bX = extrinsic(LX).b

    # the equation for u = 0
    A[1, n] = 1.0
    # new equations
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

function serial_x_in_Y(X, Y, F, tracker, cache)
    progress = cache.progress
    x0 = cache.x0
    update_x0!(x0)
    y0 = cache.y0
    y = cache.y

    LX = linear_subspace(X)
    LY = linear_subspace(Y)
    n = ambient_dim(LY)
    k = codim(LY) - codim(LX)

    dY = dim(LY)
    A, b = cache.As[dY+1], cache.bs[dY+1]

    #we loop over the points in X and check if they are contained in Y
    out = map(points(X)) do x
        update_progress!(progress, X)

        # first check
        evaluate!(y0, F, norm(x, Inf) .* x0)
        evaluate!(y, F, x)
        if norm(y, Inf) > 1e-2 * norm(y0, Inf)
            return false
        end

        # second check
        for i = 2:(k+1)
            b[i] = sum(A[i, j] * x[j] for j = 1:n)
        end
        # set up the corresponding LinearSubspace L
        E = ExtrinsicDescription(A, b; orthonormal = true)
        L = LinearSubspace(E)
        # set L as the target for homotopy continuation
        target_parameters!(tracker, L)

        U = cache.U

        # reuse U
        empty!(U)
        add!(U, x, 0)

        # add the points in Y to U after we have moved them towards L 
        for (i, p) in enumerate(points(Y))
            track!(tracker, p, 1)
            q = solution(tracker)
            _, added = add!(U, q, i)

            if !added
                return true
            end
        end

        return false
    end

    out
end
function threaded_x_in_Y(X, Y, F, tracker, cache)
    progress = cache.progress
    x0 = cache.x0
    update_x0!(x0)

    LX = linear_subspace(X)
    LY = linear_subspace(Y)
    n = ambient_dim(LY)
    k = codim(LY) - codim(LX)

    dY = dim(LY)
    A, b = cache.As[dY+1], cache.bs[dY+1]

    P = points(X)
    out = Vector{Bool}(undef, length(P))

    # Pre-allocate one tracker and buffers per thread
    nthr = Threads.nthreads()
    trackers = [deepcopy(tracker) for _ = 1:nthr]
    y0_bufs = [zeros(ComplexF64, length(cache.y0)) for _ = 1:nthr]
    y_bufs = [zeros(ComplexF64, length(cache.y)) for _ = 1:nthr]
    b_bufs = [deepcopy(b) for _ = 1:nthr]
    F_bufs = [deepcopy(F) for _ = 1:nthr]

    progress_lock = ReentrantLock()
    next_idx = Threads.Atomic{Int}(1)

    Threads.@sync begin
        for tid = 1:nthr
            let local_tracker = trackers[tid],
                local_y0 = y0_bufs[tid],
                local_y = y_bufs[tid],
                local_b = b_bufs[tid]

                local_F = F_bufs[tid]

                Threads.@spawn begin
                    while true
                        idx = Threads.atomic_add!(next_idx, 1)
                        if idx > length(P)
                            break
                        end

                        x = P[idx]

                        # first check
                        evaluate!(local_y0, local_F, norm(x, Inf) .* x0)
                        evaluate!(local_y, local_F, x)

                        result = false
                        if norm(local_y, Inf) <= 1e-2 * norm(local_y0, Inf)
                            # second check
                            for i = 2:(k+1)
                                local_b[i] = sum(A[i, j] * x[j] for j = 1:n)
                            end
                            # set up the corresponding LinearSubspace L
                            E = ExtrinsicDescription(A, local_b; orthonormal = true)
                            L = LinearSubspace(E)
                            # set L as the target for homotopy continuation
                            target_parameters!(local_tracker, L)

                            local_U = UniquePoints(x, 0)
                            add!(local_U, x, 0)

                            # add the points in Y to U after we have moved them towards L 
                            for (i, p) in enumerate(points(Y))
                                track!(local_tracker, p, 1)
                                q = solution(local_tracker)
                                _, added = add!(local_U, q, i)

                                if !added
                                    result = true
                                    break
                                end
                            end
                        end

                        lock(progress_lock) do
                            out[idx] = result
                            update_progress!(progress, X)
                        end
                    end
                end
            end
        end
    end

    return out
end

"""
    DecomposeProgress

"""

Base.@kwdef mutable struct DecomposeProgress
    progress_meter::PM.ProgressUnknown
    n_witness_sets::Int
    n_working_on::Int = 0
    degrees::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    npts::Int = 0
    step::Int = 0
end
function update_progress!(progress::DecomposeProgress)
    PM.update!(progress.progress_meter, progress.step, showvalues = showstatus(progress))
end
update_progress!(progress::Nothing, D::Vector{WitnessSet}) = nothing
update_progress_step!(progress::Nothing) = nothing
function update_progress_step!(progress::DecomposeProgress)
    progress.step += 1
    PM.update!(progress.progress_meter, progress.step, showvalues = showstatus(progress))
end
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
update_progress_n!(progress::Nothing) = nothing
function update_progress_n!(progress::DecomposeProgress)
    progress.n_working_on += 1
end
update_progress_npts!(progress::Nothing, ℓ::Int) = nothing
function update_progress_npts!(progress::DecomposeProgress, ℓ::Int) 
    progress.npts = ℓ
end
finish_progress!(progress::DecomposeProgress) =
    PM.finish!(progress.progress_meter, showvalues = showstatus(progress))
function showstatus(progress::DecomposeProgress)
    ℓ = progress.npts
    if ℓ > 0
        text = [("Decomposing witness set", "$(progress.n_working_on)/$(progress.n_witness_sets) ($ℓ points left to classify)")]
    else
        text = [("Decomposing witness set", "$(progress.n_working_on)/$(progress.n_witness_sets)")]
    end
    degs = progress.degrees
    k = sort(collect(keys(degs)), rev = true)
    for key in k
        push!(
            text,
            ("Degrees of components of dim. $key", "$(join(map(string, degs[key]), ','))"),
        )
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
    threading::Bool = true,
)


    P = points(W)
    L = linear_subspace(W)
    G = system(W)
    n = ambient_dim(L)

    decomposition = Vector{WitnessSet}()

    if dim(L) < n

        MS = MonodromySolver(G, L; compile = false, options = options)

        res = monodromy_solve(
            MS,
            P,
            L,
            seed;
            threading = threading,
            show_progress = show_monodromy_progress,
        )

        if warning && (trace(res) > options.trace_test_tol)
            @warn "Trying to decompose non-complete set of witness points for codimension $(dim(L)) (trace test failed). Will try to compute the missing points. The output will contain all components, for which the trace test was successfull."
        end

        iter = 0
        non_complete_points = solutions(res)
        d = length(non_complete_points) # the total degree (i.e., number of points on all irreducible components)
        non_complete_orbits = Vector{Set{Int}}()

        while !isempty(non_complete_points)
            update_progress!(progress)
            update_progress_step!(progress)

            iter += 1
            if iter > max_iters
                break
            end

            # for safety an additional monodromy
            res = monodromy_solve(
                MS,
                non_complete_points,
                L,
                seed;
                threading = threading,
                show_progress = show_monodromy_progress,
            )
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
                update_progress!(progress)
                update_progress_step!(progress)

                P_orbit = non_complete_points[collect(orbit)]
                res_orbit = monodromy_solve(
                    MS,
                    P_orbit,
                    L,
                    seed;
                    threading = threading,
                    show_progress = show_monodromy_progress,
                )

                if trace(res_orbit) < options.trace_test_tol

                    # We do not want to add orbits of length 1 in the beginning. Even if they are on an irreducible component of degree > 1, they tend to have small trace.
                    if length(orbit) > 1 || iter ≥ 5
                        W_new = WitnessSet(G, L, P_orbit)

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

        update_progress_npts!(progress, 0)
    else
        for p in P
            push!(decomposition, WitnessSet(G, L, [p]))
        end
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
* `show_progress = true`: indicate whether the progress of the computation should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress of the monodromy computation should be displayed.
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_iters = 50`: maximal number of iterations for the decomposition step.
* `warning = true`: if `true` prints a warning when the [`trace_test`](@ref) fails. 
* `threading = true`: enables multiple threads. 
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
    Ws::Union{Vector{WitnessSet{T1,T2,Vector{T3}}},Vector{WitnessSet}};
    show_progress::Bool = true,
    show_monodromy_progress::Bool = false,
    monodromy_options::MonodromyOptions = MonodromyOptions(; trace_test_tol = 1e-10),
    max_iters::Int = 50,
    warning::Bool = true,
    threading::Bool = true,
    seed = nothing,
) where {T1,T2,T3<:Number}

    if isnothing(seed)
        seed = rand(UInt32)
    end
    Random.seed!(seed)

    options = decompose_with_monodromy_options(monodromy_options)
    out = Vector{WitnessSet}()

    if isempty(Ws)
        return out
    end

    c = length(Ws)
    n = ambient_dim(linear_subspace(Ws[1]))

    if show_progress
        progress = DecomposeProgress(
            progress_meter =
            PM.ProgressUnknown(
                dt = 0.1,
                desc = "Decomposing $c witness sets",
                enabled = true,
                spinner = true,
            ),
            n_witness_sets = c
        )
    else
        progress = nothing
    end

    for W in Ws
        if ModelKit.degree(W) > 0

            update_progress_n!(progress)
            update_progress_step!(progress)

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
    header = "\n Numerical irreducible decomposition with $total components"
    println(io, header)

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
        header = headers,
        tf = PrettyTables.tf_unicode_rounded,
        alignment = :c,
        header_crayon = PrettyTables.Crayon(bold = false),
        border_crayon = PrettyTables.Crayon(faint = true),
    )
end

"""
    numerical_irreducible_decomposition(F::System; options...)

Computes the numerical irreducible of the variety defined by ``F=0``. 

### Options

* `show_progress = true`: indicate whether the progress of the computation should be displayed.
* `show_monodromy_progress = false`: indicate whether the progress of the monodromy computation should be displayed.
* `sorted = true`: the polynomials in F will be sorted by degree in decreasing order. 
* `max_codim`: the maximal codimension until which witness supersets should be computed.
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `monodromy_options`: [`MonodromyOptions`](@ref) for [`monodromy_solve`](@ref).
* `max_iters = 50`: maximal number of iterations for the decomposition step.
* `warning = true`: if `true` prints a warning when the [`trace_test`](@ref) fails. 
* `threading = true`: enables multiple threads. 
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
"""
function numerical_irreducible_decomposition(
    F::System;
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(;
        max_endgame_steps = 100,
        max_endgame_extended_steps = 100,
        sing_accuracy = 1e-10,
    ),
    show_monodromy_progress::Bool = false,
    monodromy_options::MonodromyOptions = MonodromyOptions(; trace_test_tol = 1e-10),
    max_iters::Int = 50,
    sorted::Bool = true,
    max_codim::Union{Int,Nothing} = nothing,
    warning::Bool = true,
    threading::Bool = true,
    seed = nothing,
    kwargs...,
)

    !isnothing(seed) && Random.seed!(seed)

    Ws = regeneration(
        F;
        sorted = sorted,
        max_codim = max_codim,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
        threading = threading,
        seed = nothing,
        kwargs...,
    )
    if isnothing(Ws)
        return nothing
    end
    dec = decompose(
        Ws;
        monodromy_options = monodromy_options,
        max_iters = max_iters,
        show_monodromy_progress = show_monodromy_progress,
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

Numerical irreducible decomposition with 1 components
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
