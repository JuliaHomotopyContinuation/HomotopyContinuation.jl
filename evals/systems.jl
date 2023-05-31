using HomotopyContinuation

Base.@kwdef struct EvalTestSystem
    system::System
    nsolutions::Int
end

function lines_clebsch_cubic()
    @var x y z t a[1:3] b[1:3]
    f =
        81 * (x^3 + y^3 + z^3) -
        189 * (x^2 * y + x^2 * z + x * y^2 + x * z^2 + y^2 * z + y * z^2) +
        54 * x * y * z +
        126 * (x * y + x * z + y * z) - 9 * (x^2 + y^2 + z^2) - 9 * (x + y + z) + 1
    fab = subs(f, [x; y; z] => a + t * b)
    _, C = exponents_coefficients(fab, [t])
    F = subs(
        C,
        [a[3]; b[3]] => [-(7 + a[1] + 3 * a[2]) / 5; -(11 + 3 * b[1] + 5 * b[2]) / 7],
    )

    system = System(F)
    nsolutions = 27
    return EvalTestSystem(system, nsolutions)
end



function simple_solve(tracker, starts; threaded = false)
    if threaded
        S = collect(enumerate(starts))
        N = length(S)
        results = Vector{PathResult}(undef, N)
        trackers = [tracker; [deepcopy(tracker) for _ = 2:Threads.nthreads()]]
        Threads.@threads for (k, start) in S
            results[k] = track(trackers[Threads.threadid()], start)
        end
        return results
    else
        return map(x -> track(tracker, x), start)
    end
end

function run_single_system_benchmark(
    evaltestsystem;
    threaded = false,
    gamma = cis(2π * rand()),
)
    F = evaltestsystem.system
    tracker, start = total_degree(F; gamma = gamma)


    timedresult = @timed simple_solve(tracker, start; threaded = threaded)
    result = timedresult.value
    found_solutions = count(r -> is_success(r), result)
    correct_nsolutions = found_solutions == evaltestsystem.nsolutions
    time = timedresult.time

    return Dict(
        :found_solutions => found_solutions,
        :correct_nsolutions => correct_nsolutions,
        :time => time,
        :result => result,
    )
end


function run_system_benchmark(evaltestsystem; nruns = 10)
    gammas = [cis(2π * rand()) for _ = 1:nruns]
    bench_results =
        map(gamma -> run_single_system_benchmark(evaltestsystem; gamma = gamma), gammas)

    percentage_correct = count(r -> r[:correct_nsolutions], bench_result) / nruns

    return Dict(:percentage_correct => percentage_correct, :bench_results => bench_results)
end



sys1 = lines_clebsch_cubic()

run_system_benchmark(sys1)

F = lines_clebsch_cubic().system

tracker, start = total_degree(F, tracker_options = (max_steps = 200,))

S = collect(start)

@time simple_solve(tracker, S; threaded = true);


s = first(start)

@timed track(tracker, S[4])



R = map(x -> track(tracker, x), start)




sort(map(r -> r.accepted_steps, R))
