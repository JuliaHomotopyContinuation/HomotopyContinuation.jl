"""
    OverdeterminedTracker

An `OverdeterminedTracker` allows to solve overdetermined systems by
first squaring them up, then solving them with a `PolyhedralTracker` or `PathTracker`
and finally to find all actual solutions.
"""
struct OverdeterminedTracker{
    T<:AbstractPathTracker,
    AV<:AbstractVector,
    S<:HomotopyWithCache,
    NC<:NewtonCorrectorCache,
} <: AbstractPathTracker
    tracker::T
    # stuff to verify containment
    system::S
    y::AV
    jacobian::JacobianMonitor{Float64}
    newton::NC
end

function construct_tracker(prob::OverdeterminedProblem, start_solutions; kwargs...)
    tracker = construct_tracker(prob.problem, start_solutions)
    y = copy(path_tracker_state(tracker).solution)
    system = HomotopyWithCache(ConstantHomotopy(prob.target_system), y, 0.0)
    J = JacobianMonitor(jacobian(system, y, 0.0))
    newton_corrector = NewtonCorrectorCache(system, y, 0.0)

    OverdeterminedTracker(tracker, system, y, J, newton_corrector)
end


seed(OT::OverdeterminedTracker) = seed(OT.tracker)
function PathResult(OT::OverdeterminedTracker, x, path_number = nothing; kwargs...)
    PathResult(OT.tracker, x, path_number; kwargs...)
end

result_type(OT::OverdeterminedTracker) = result_type(OT.tracker)
prepare!(OT::OverdeterminedTracker, S) = prepare!(OT.tracker, S)
function track!(
    OT::OverdeterminedTracker,
    x;
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
)
    retcode = track!(
        OT.tracker,
        x;
        accuracy = accuracy,
        max_corrector_iters = max_corrector_iters,
    )
    is_success(retcode) || return retcode

    # We take the solution and try an overdetermined Newton method to see whether this also
    # converges
    s = solution(OT.tracker)
    result = newton!(
        OT.y,
        OT.system,
        s,
        0.0,
        OT.jacobian,
        norm(OT.tracker),
        OT.newton;
        tol = min_accuracy(OT.tracker),
        # Make sure to evaluate the Jacobian only *once*. Otherwise it can happen at singuar
        # solutions that we bounce away from a good solution.
        full_steps = 1,
        max_iters = 2,
        double_64_evaluation = false,
    )

    # overwrite the solution information used for constructing a `PathResult`
    state = path_tracker_state(OT)
    state.solution_cond = cond!(OT.jacobian, InfNorm())
    original_residual = state.solution_residual
    state.solution_residual = LA.norm(OT.newton.r, InfNorm())
    state.solution_accuracy = result.accuracy

    if is_converged(result)
        state.solution .= OT.y
        state.status = PathTrackerStatus.success
    else
        if state.solution_cond > 1e8
            evaluate!(OT.newton.r, OT.system, s, 0.0)
            state.solution_residual = LA.norm(OT.newton.r, InfNorm())
            if state.solution_residual / original_residual < 1e5
                state.status = PathTrackerStatus.success
            else
                state.status = PathTrackerStatus.excess_solution
            end
        else
            state.status = PathTrackerStatus.excess_solution
        end
    end
end

function track(
    tracker::OverdeterminedTracker,
    x;
    path_number::Union{Int,Nothing} = nothing,
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
    details::Symbol = :default,
)
    track!(tracker, x; accuracy = accuracy, max_corrector_iters = max_corrector_iters)
    PathResult(tracker, x, path_number; details = details)
end

path_tracker_state(OT::OverdeterminedTracker) = path_tracker_state(OT.tracker)
path_tracker_state(PT::PolyhedralTracker) = path_tracker_state(PT.generic_tracker)
path_tracker_state(T::PathTracker) = T.state
