export WitnessSet,
    witness_set, linear_subspace, system, dim, codim, trace_test, is_irreducible, membership

"""
    WitnessSet(F, L, S)

Store solutions `S` of the polynomial system `F(x) = L(x) = 0` into a witness set.
"""
struct WitnessSet{
    S<:AbstractSystem,
    Sub<:AbstractSubspace,
    R<:Union{Vector{ComplexF64},PathResult},
}
    F::S
    L::Sub
    # only non-singular results or solutions
    R::Vector{R}
    projective::Bool
    # is_irreducible is being set by decompose
    is_irreducible::Union{Nothing,Bool}
end

function WitnessSet(
    F::AbstractSystem,
    L::LinearSubspace,
    R;
    projective::Bool = is_linear(L) && is_homogeneous(System(F)),
    is_irreducible::Union{Nothing,Bool} = nothing,
)
    WitnessSet(F, L, R, projective, is_irreducible)
end

"""
    system(W::WitnessSet)

Get the system stored in `W`.
"""
system(W::WitnessSet) = W.F

"""
    linear_subspace(W::WitnessSet)

Get the linear subspace stored in `W`.
"""
linear_subspace(W::WitnessSet) = W.L

"""
    solutions(W::WitnessSet)

Get the solutions stored in `W`.
"""
solutions(W::WitnessSet{A,B,PathResult}) where {A,B} = solutions(W.R)
solutions(W::WitnessSet{A,B,Vector{ComplexF64}}) where {A,B} = W.R

points(W::WitnessSet) = solutions(W)

"""
    results(W::WitnessSet)

Get the results stored in `W`.
"""
results(W::WitnessSet{<:Any,<:Any,PathResult}) = W.R

"""
    degree(W::WitnessSet)

Returns the degree of the witness set `W`. This equals the number of solutions stored.
"""
ModelKit.degree(W::WitnessSet) = length(W.R)

"""
    dim(W::WitnessSet)

The dimension of the algebraic set encoded by the witness set.
"""
dim(W::WitnessSet) = codim(W.L)

"""
    codim(W::WitnessSet)

The dimension of the algebraic set encoded by the witness set.
"""
codim(W::WitnessSet) = dim(W.L)

function Base.show(io::IO, W::WitnessSet)
    print(io, "Witness set for dimension $(dim(W)) of degree $(degree(W))")
end

copy(W::WitnessSet) =  WitnessSet(copy(W.F), copy(W.L), copy(W.R), copy(W.projective), copy(W.is_irreducible))



"""
    is_irreducible(W::WitnessSet)

Return `true`, if it was computed that `W` is irreducible, `false` if it was computed that `W` is not irreducible, and `:undecided` otherwise.

### Example 
```julia
julia> @var x[1:2]
julia> f = System([sum(x.^2) - 1])
julia> W = witness_set(f; dim = 1);
julia> is_irreducible(W)
:undecided

julia> dec = decompose(W);
julia> W_irred = first(dec);
julia> is_irreducible(W_irred)
true
end
```
"""
function is_irreducible(W::WitnessSet)
    if isnothing(W.is_irreducible)
        :undecided
    else
        W.is_irreducible
    end
end


### Construct witness sets
"""
    witness_set(F; codim = nvariables(F) - length(F), dim = nothing, options...)

Compute a [`WitnessSet`](@ref) for `F` in the given dimension (resp. codimension)
by sampling a random (affine) linear subspace.
After constructing the system this calls [`solve`](@ref) with the provided `options`.

    witness_set(F, L; options...)

Compute [`WitnessSet`](@ref) for `F` and the (affine) linear subspace `L`.

    witness_set(W::WitnessSet, L; options...)

Compute a new [`WitnessSet`](@ref) with the (affine) linear subspace `L` by moving
the linear subspace stored in `W` to `L`.

### Example
```julia-repl
julia> @var x y;
julia> F = System([x^2 + y^2 - 5], [x, y])
System of length 1
 2 variables: x, y

 -5 + x^2 + y^2

julia> W = witness_set(F)
Witness set for dimension 1 of degree 2
```
"""
witness_set(F::System, args...; compile = COMPILE_DEFAULT[], kwargs...) =
    witness_set(fixed(F; compile = compile), args...; kwargs...)
function witness_set(
    F::AbstractSystem;
    target_parameters = nothing,
    dim = nothing,
    codim = nothing,
    options...,
)
    f = System(F)
    projective = is_homogeneous(f)
    if isnothing(dim) && isnothing(codim)
        dim = corank(F) - projective
    elseif !isnothing(codim) && projective
        codim += 1
    end
    L = rand_subspace(size(F, 2); dim = codim, codim = dim, affine = !projective)
    witness_set(
        F,
        L;
        target_parameters = target_parameters,
        projective = projective,
        options...,
    )
end

function witness_set(F::AbstractSystem, L::LinearSubspace; projective = nothing, options...)
    res = solve(F; target_subspace = L, options...)
    if isnothing(projective)
        WitnessSet(F, L, results(res; only_nonsingular = true))
    else
        WitnessSet(F, L, results(res; only_nonsingular = true); projective = projective)
    end
end

corank(F::AbstractSystem; kwargs...) = nvariables(F) - LA.rank(F)
function LinearAlgebra.rank(F::AbstractSystem; target_parameters = nothing)
    m, n = size(F)
    u = zeros(ComplexF64, m)
    U = zeros(ComplexF64, m, n)
    x = randn(ComplexF64, n)
    evaluate_and_jacobian!(u, U, F, x, target_parameters)
    LA.rank(U)
end


### Move witness sets around
function witness_set(W::WitnessSet, L::LinearSubspace; options...)
    if W.projective && !is_linear(L)
        error(
            "The given space is an affine linear subspace (``b ≠ 0``). " *
            " Expected a linear subspace since the given witness set is projective.",
        )
    end
    res = solve(W.F, W.R; start_subspace = W.L, target_subspace = L, options...)
    WitnessSet(
        W.F,
        L,
        results(res; only_nonsingular = true);
        projective = is_linear(L) && W.projective,
    )
end

### Membership
Base.@kwdef mutable struct MembershipProgress
    progress_meter::PM.Progress
    current_task::Int = 0
    ntasks::Int = 0
end
MembershipProgress(progress_meter::PM.Progress) = MembershipProgress(progress_meter = progress_meter)
update_progress_tasks!(progress::Nothing, i::Int, m::Int) = nothing
function update_progress_tasks!(progress::MembershipProgress, i::Int, m::Int)
    progress.current_task = i
    progress.ntasks = m
    PM.update!(progress.progress_meter, i; showvalues = [("Points checked", "$(progress.current_task)/$(progress.ntasks)")])
end
update_progress!(progress::Union{Nothing,MembershipProgress}, W::Nothing) = nothing

"""
    MembershipCache

A cache for [`membership`](@ref).
"""
mutable struct MembershipCache
    A::Matrix
    b::Vector
    x0::Vector
    y0::Vector
    y::Vector

    endgame_options::EndgameOptions
    tracker_options::TrackerOptions

    progress::Union{MembershipProgress,Nothing}
end

function MembershipCache(W, EO, TO, progress)
    F = system(W)
    m, n = size(F)
    i = codim(W)
    A0 = zeros(ComplexF64, n - i, n)
    A = LA.svd(A0).Vt # need to orthonormalize A
    b = zeros(ComplexF64, n - i)
    x0 = zeros(ComplexF64, n)
    y0 = zeros(ComplexF64, m)
    y = zeros(ComplexF64, m)

    MembershipCache(A, b, x0, y0, y, EO, TO, progress)
end
function update_x0!(x0)
    for i = 1:length(x0)
        x0[i] = randn(ComplexF64)
    end
    LA.normalize!(x0)
    x0
end

"""
    membership(P::Vector{Vector{T}}, W::WitnessSet) where {T <: Number}

Returns a boolean vector indicating whether the points of P are contained in X.

### Options

* `show_progress = true`: indicate whether a progress bar should be displayed.
* `tracker_options`: [`TrackerOptions`](@ref) for the [`Tracker`](@ref).
* `endgame_options`: [`EndgameOptions`](@ref) for the [`EndgameTracker`](@ref).
* `atol = 1e-14` and `rtol = sqrt(eps())`: a point `y` is considered equal to `x` when the distance between `x`and `y` is smaller than `max(atol, norm(x, Inf) * rtol).`
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)
"""
function membership(P::Vector{Vector{T}}, W::WitnessSet; 
                    show_progress::Bool = true,
                    tracker_options = TrackerOptions(),
                    endgame_options = EndgameOptions(;
                        max_endgame_steps = 100,
                        max_endgame_extended_steps = 100,
                        sing_cond = 1e12,
                    ),
                    kwargs...) where {T <: Number}

    # progress bar
    if show_progress
        desc = "Membership test"
        barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
        progress_meter = ProgressMeter.Progress(length(P); dt = 0.2, desc = desc, barlen = barlen, output = stdout)
        progress = MembershipProgress(progress_meter)
    else
        progress = nothing
    end

    cache = MembershipCache(W, endgame_options, tracker_options, progress)
    is_contained(P, W, system(W), cache; kwargs...)
end
function membership(p::Vector{T}, W::WitnessSet; kwargs...) where {T <: Number}
    out = membership([p], W)
    first(out)
end
function is_contained(P::Vector{Vector{T}}, Y::WitnessSet, F, cache; threading::Bool = Threads.nthreads() > 1, kwargs...) where {T <: Number}
    
    # set up homotopy
    tracker_options = cache.tracker_options
    endgame_options = cache.endgame_options
    LY = linear_subspace(Y)
    Hom = linear_subspace_homotopy(F, LY, LY)
        tracker = EndgameTracker(
            Hom;
            tracker_options = tracker_options,
            options = endgame_options,
        )

    # tracking
    if threading
        out = threaded_x_in_Y(P, Y, F, tracker, cache; kwargs...)
    else
        out = serial_x_in_Y(P, Y, F, tracker, cache; kwargs...)
    end

    out
end
function serial_x_in_Y(P, Y, F, tracker, cache; atol = 1e-14, rtol = sqrt(eps()))

    progress = cache.progress
    x0 = cache.x0
    update_x0!(x0)
    y0 = cache.y0
    y = cache.y
    l_X = length(P)

    if first(cache.b) isa Number
        A, b =  cache.A, cache.b
    else
        LY = linear_subspace(Y)
        dY = dim(LY)
        A, b = cache.A[dY+1], cache.b[dY+1]
    end
   
    # Pre-allocate output
    out = Vector{Bool}(undef, l_X)
    idx = 0

    #we loop over the points in X and check if they are contained in Y
    out = map(P) do x
        idx += 1
        update_progress_tasks!(progress, idx, l_X)

        # first check
        evaluate!(y0, F, norm(x, Inf) .* x0)
        evaluate!(y, F, x)
        if norm(y, Inf) > 1e-2 * norm(y0, Inf)
            return false
        end

        # second check
        LA.mul!(b, A, x)
        # set up the corresponding LinearSubspace L
        E = ExtrinsicDescription(A, b; orthonormal = true)
        L = LinearSubspace(E)
        # set L as the target for homotopy continuation
        target_parameters!(tracker, L)

        rad = max(atol, norm(x, Inf) * rtol)

        # add the points in Y to U after we have moved them towards L 
        for p in points(Y)
            track!(tracker, p, 1)
            q = solution(tracker)
            d = distance(q, x, InfNorm())
            if  d < rad
                return true
            end
        end

        return false
    end

    out
end
function threaded_x_in_Y(P, Y, F, tracker, cache; atol = 1e-14, rtol = sqrt(eps()))

    progress = cache.progress
    x0 = cache.x0
    update_x0!(x0)
    y0 = cache.y0
    y = cache.y
    l_X = length(P)

    if first(cache.b) isa Number
        A, b =  cache.A, cache.b
    else
        LY = linear_subspace(Y)
        dY = dim(LY)
        A, b = cache.A[dY+1], cache.b[dY+1]
    end

    # Pre-allocate output
    out = Vector{Bool}(undef, l_X)

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
                        if idx > l_X
                            break
                        end

                        x = P[idx]
                        lock(progress_lock) do
                            update_progress_tasks!(progress, idx, l_X)
                        end

                        # first check
                        evaluate!(local_y0, local_F, norm(x, Inf) .* x0)
                        evaluate!(local_y, local_F, x)

                        result = false
                        if norm(local_y, Inf) <= 1e-2 * norm(local_y0, Inf)
                            # second check
                            LA.mul!(local_b, A, x)
                            # set up the corresponding LinearSubspace L
                            E = ExtrinsicDescription(A, local_b; orthonormal = true)
                            L = LinearSubspace(E)
                            # set L as the target for homotopy continuation
                            target_parameters!(local_tracker, L)

                            rad = max(atol, norm(x, Inf) * rtol)

                            # add the points in Y to U after we have moved them towards L 
                            for p in points(Y)
                                track!(local_tracker, p, 1)
                                q = solution(local_tracker)
                                if distance(q, x, InfNorm()) < rad
                                    result = true
                                    break
                                end
                            end
                        end

                        lock(progress_lock) do
                            out[idx] = result
                        end
                    end
                end
            end
        end
    end

    return out
end

"""
    trace_test(W::WitnessSet; options...)

Performs a trace test [^LRS18] to verify whether the given witness set `W` is complete.
Returns the trace of the witness set which should be theoretically be 0 if `W` is complete.
Due to floating point arithmetic this is not the case, thus is has to be manually checked
that the trace is sufficiently small.
Returns `nothing` if the trace test failed due to path tracking failures.
The `options` are the same as for calls to [`witness_set`](@ref).

```julia-repl
julia> @var x y;
julia> F = System([x^2 + y^2 - 5], [x, y])
System of length 1
 2 variables: x, y

 -5 + x^2 + y^2

julia> W = witness_set(F)
Witness set for dimension 1 of degree 2

julia> trace = trace_test(W)
9.981960497718987e-16
```

[^LRS18]: Leykin, Anton, Jose Israel Rodriguez, and Frank Sottile. "Trace test." Arnold Mathematical Journal 4.1 (2018): 113-125.
APA

"""
function trace_test(W₀::WitnessSet; options...)
    L₀ = linear_subspace(W₀)
    F = system(W₀)
    S₀ = solutions(W₀)
    # if we are in the projective setting, we need to make sure that
    # all solutions are on the same affine chart
    # Therefore make the affine chart now
    if W₀.projective
        F = on_affine_chart(F)
        s₀ = sum(s -> set_solution!(s, F, s), S₀)
    else
        s₀ = sum(S₀)
    end

    v = randn(ComplexF64, codim(L₀))
    L₁ = translate(L₀, v)
    L₋₁ = translate(L₀, -v)

    R₁ = solve(F, S₀; start_subspace = L₀, target_subspace = L₁, options...)
    nsolutions(R₁) == degree(W₀) || return nothing

    R₋₁ = solve(F, S₀; start_subspace = L₀, target_subspace = L₋₁, options...)
    nsolutions(R₋₁) == degree(W₀) || return nothing

    s₁ = sum(solutions(R₁))
    s₋₁ = sum(solutions(R₋₁))

    M = [s₋₁ s₀ s₁; 1 1 1]
    singvals = LA.svdvals(M)
    trace = singvals[3] / singvals[1]

    trace
end

