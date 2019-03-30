export new_solve

function new_solve(args...; threading=true, report_progress=true, kwargs...)
    tracker, start_solutions = pathtracker_startsolutions(args...; kwargs...)
    solve(tracker, start_solutions; threading=threading, report_progress=report_progress)
end

function solve(tracker::PathTracker, start_solutions; kwargs...)
    results = Vector{result_type(tracker)}(undef, length(start_solutions))
    track_paths!(results, tracker, start_solutions; kwargs...)
    results
end

function track_paths!(results, tracker, start_solutions; threading=true, report_progress=true, details_level=1)
    track_paths!(results, tracker, start_solutions, threading, report_progress, details_level)
end
function track_paths!(results, tracker, start_solutions, threading, report_progress, details_level)
    n = length(results)

    if report_progress
        progress = ProgressMeter.Progress(n, 0.1, "Tracking $n paths... ")
    else
        progress = nothing
    end

    nthreads = Threads.nthreads()
    if threading && nthreads > 1
        # TODO: We can probably also do this better, but for now we have to collect
        # to support indexing
        S = collect(start_solutions)

        batch_size = 32 * nthreads
        ranges = partition_work(1:min(batch_size, n), nthreads)
        trackers = Threads.resize_nthreads!([tracker])
        batch_tracker = BatchTracker(results, trackers, ranges, S, details_level)

        k = 1
        while k ≤ n
            partition_work!(batch_tracker.ranges, k:min(k+batch_size-1, n), nthreads)
            ccall(:jl_threading_run, Ref{Cvoid}, (Any,), batch_tracker)
            k += batch_size
            update_progress!(progress, results, min(k - 1, n))
        end
    else
        for (k, s) in enumerate(start_solutions)
            results[k] = track(tracker, s, 1.0, k; details_level=details_level)
            k % 16 == 0 && update_progress!(progress, results, k)
        end
    end
    results
end

function update_progress!(progress, results, N)
    ProgressMeter.update!(progress, N, showvalues=((:tracked, N),))
    nothing
end
update_progress!(::Nothing, results, N) = nothing

mutable struct BatchTracker{Tracker<:PathTracker, V, R} <: Function
    results::Vector{R}
    trackers::Vector{Tracker}
    ranges::Vector{UnitRange{Int}}
    start_solutions::V
    details_level::Int
end

function (batch::BatchTracker)()
    tid = Threads.threadid()
    track_batch!(batch.results, batch.trackers[tid],
                 batch.ranges[tid], batch.start_solutions, batch.details_level)
end
function track_batch!(results, pathtracker, range, starts, details_level)
    for k in range
        results[k] = track(pathtracker, starts[k], 1.0, k; details_level=details_level)
    end
    results
end
