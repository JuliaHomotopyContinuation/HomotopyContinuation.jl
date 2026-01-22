export IntrinsicSubspaceTracker

"""
    IntrinsicSubspaceTracker <: AbstractPathTracker

"""
struct IntrinsicSubspaceTracker{H1<:IntrinsicSubspaceStiefelHomotopy, H2<:IntrinsicSubspaceOffsetHomotopy, M} <: AbstractPathTracker
    stiefel_tracker::Tracker{H1,M}
    offset_tracker::EndgameTracker{H2,M}
end
IntrinsicSubspaceTracker(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
) = IntrinsicSubspaceTracker(fixed(F; compile = compile), start, target; kwargs...)
function IntrinsicSubspaceTracker(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
)
    IntrinsicSubspaceTracker(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
        kwargs...
    )
end

function IntrinsicSubspaceTracker(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64};
    endgame_options = EndgameOptions(),
    tracker_options = TrackerOptions(),
    affine::Bool = true,
    kwargs...
)
    unsupported_kwargs(kwargs)
    if affine 
        H1 = IntrinsicSubspaceStiefelHomotopy(system, start, target; affine = true)
        H2 = IntrinsicSubspaceOffsetHomotopy(system, start, target)
        T1 = Tracker(H1; options = tracker_options)
        T2 = Tracker(H2; options = tracker_options)
        return IntrinsicSubspaceTracker(T1, EndgameTracker(T2; options = endgame_options))
    else
        H = IntrinsicSubspaceStiefelHomotopy(system, start, target; affine = false)
        return EndgameTracker(Tracker(H))
    end
end

function track(
    tracker::IntrinsicSubspaceTracker,
    start_solution::AbstractVector{ComplexF64};
    path_number::Union{Nothing,Int} = nothing,
    debug::Bool = false,
)
   
    # The tracker works in two stages
    # first it runs a Stiefel homotopy (Ax + a) -> (Bx + a)
    # then it runs an offset homotopy (Bx + a) -> (Bx + b)
    retcode = track!(
            tracker.stiefel_tracker,
            start_solution;
            debug = debug,
        )


    @unpack μ, ω = tracker.stiefel_tracker.state
    state = tracker.stiefel_tracker.state

    if !is_success(retcode)
        t = real(state.t)
        x = get_solution(tracker.stiefel_tracker.homotopy, state.x, t) 
        return PathResult(
            return_code = :stiefel_homotopy_failed,
            solution = x,
            start_solution = start_solution,
            t = t,
            accuracy = state.accuracy,
            singular = false,
            condition_jacobian = NaN,
            residual = NaN,
            winding_number = nothing,
            last_path_point = (x, t),
            valuation = nothing,
            ω = state.ω,
            μ = state.μ,
            path_number = path_number,
            extended_precision = state.extended_prec,
            accepted_steps = state.accepted_steps,
            rejected_steps = state.rejected_steps,
            extended_precision_used = state.used_extended_prec,
        )
    end

    x = get_solution(tracker.stiefel_tracker.homotopy, state.x, 0.0) 
    r = track(
        tracker.offset_tracker,
        x;
        path_number = path_number,
        debug = debug,
    )
    r.start_solution = start_solution
    # report accurate steps
    r.accepted_steps += tracker.stiefel_tracker.state.accepted_steps
    r.rejected_steps += tracker.stiefel_tracker.state.rejected_steps

    r
end