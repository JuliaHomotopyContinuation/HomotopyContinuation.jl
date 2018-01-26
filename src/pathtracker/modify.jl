export setup_pathtracker!

"""
    setup_pathtracker!(tracker, x0, s_start, s_end)

Reset the given pathtracker `tracker` and set it up to track `x0` form `s_start` to
`s_end`.
"""
function setup_pathtracker!(tracker::Pathtracker{Low}, x0, s_start, s_end) where Low
    # TODO: we should carry over the original options
    setprecisionvalues!(tracker, x0)
    tracker.startvalue .= tracker.low.x
    #tracker.startvalue = copy(tracker.low.x)

    tracker.iter = 0
    tracker.t = 1.0
    tracker.steplength = tracker.options.initial_steplength
    tracker.s = convert(Complex{Low}, s_start)
    tracker.sdiff = convert(Complex{Low}, s_end) - tracker.s
    tracker.ds = tracker.steplength * tracker.sdiff
    tracker.snext = tracker.s
    tracker.step_sucessfull = false
    tracker.consecutive_successfull_steps = 0
    tracker.status = :default

    tracker
end

function setprecisionvalues!(tracker::Pathtracker{Low, High}, x0::AbstractVector) where {Low, High}
    N = length(tracker.low.x)
    n = length(x0)
    if n == N - 1
        tracker.low.x[2:end] .= x0
        tracker.low.x[1] = one(Low)
    elseif n == N
        tracker.low.x .= x0
    else
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    tracker.low.xnext .= tracker.low.x
    tracker.high.x .= tracker.low.x #this automatically converts to High
    tracker.high.xnext .= tracker.high.x
    tracker.usehigh = false
end

function setprecisionvalues!(tracker::Pathtracker{Low, High}, x0::AbstractVector{Complex{High}}) where {Low, High}
    N = length(tracker.high.x)
    if length(x0) == N - 1
        tracker.high.x[2:end] .= x0
        tracker.high.x[1] = one(High)
    elseif length(x0) == N
        tracker.high.x .= x0
    else
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    tracker.high.xnext .= tracker.high.x
    tracker.low.x .= tracker.high.x #this automatically converts to Low
    tracker.low.xnext .= tracker.low.x

    tracker.usehigh = true
end
