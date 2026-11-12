# Run the same script with --project pointing at each checkout being compared.
# Reports complete, warmed tracking runs; it does not benchmark certification.
using HomotopyContinuation, LinearAlgebra, Random, Printf
const HC = HomotopyContinuation

function tracking_cases()
    @var x y p q
    coupled = (
        name = "coupled_roots",
        F = System([x^2 - p, y^2 - x]; variables = [x, y], parameters = [p]),
        p = [1.0],
        q = [4.0],
        starts = [[1, 1], [1, -1], [-1, im], [-1, -im]],
        targets = [[2, sqrt(2)], [2, -sqrt(2)], [-2, im * sqrt(2)], [-2, -im * sqrt(2)]],
    )
    @var z[1:8] a[1:8]
    M = Matrix{Float64}(I, 8, 8) + 0.1randn(MersenneTwister(42), 8, 8)
    target = collect(range(1.1, 1.8; length = 8))
    dense = (
        name = "dense_cubics",
        F = System(M * (z .^ 3 - a); variables = z, parameters = a),
        p = ones(8),
        q = target,
        starts = [ones(8)],
        targets = [cbrt.(target)],
    )
    extended = (
        name = "rounded_jacobian",
        F = System(
            [x + y - 1, x + y - 1 + p * ((1 + q) * y - 1)];
            variables = [x, y],
            parameters = [p, q],
        ),
        p = [2.0^-60, 1],
        q = [2.0^-60, 0],
        starts = [[0.5, 0.5]],
        targets = [[0.0, 1.0]],
    )
    [coupled, dense, extended, steiner_case()]
end

function steiner_case()
    # Reuse the existing difficult-path fixture, without running its test suite
    # or maintaining a second copy of the numerical data and polynomial system.
    source = read(joinpath(@__DIR__, "../test/test_cases/steiner_higher_prec.jl"), String)
    marker = "@testset \"Steiner High Precision\" begin"
    first_data = nextind(source, last(findfirst(marker, source)))
    last_data = prevind(source, first(findfirst("    @testset \"Steiner - AD:", source)))
    body = source[first_data:last_data]
    fixture = include_string(
        @__MODULE__,
        "function steiner_fixture()\n" *
        body *
        "\n(name=\"steiner_precision\", F=F, p=p, q=q, starts=[s_p], targets=[s_q])\nend",
    )
    Base.invokelatest(fixture)
end

function batch!(tracker, starts, repetitions)
    for _ = 1:repetitions, start in starts
        track!(tracker, start)
    end
end

function benchmark_case(case, compile; samples = 5, repetitions = 10)
    H = ParameterHomotopy(case.F, case.p, case.q; compile)
    tracker = Tracker(H)
    batch!(tracker, case.starts, 2)
    times = Float64[]
    for _ = 1:samples
        before = time_ns()
        batch!(tracker, case.starts, repetitions)
        push!(times, (time_ns() - before) * 1e-9 / (repetitions * length(case.starts)))
    end
    successes = accepted = rejected = 0
    error = residual = 0.0
    u = zeros(HC.DoubleDouble.ComplexDF64, size(H, 1))
    for (start, target) in zip(case.starts, case.targets)
        result = track(tracker, start)
        successes += is_success(result)
        accepted += accepted_steps(result)
        rejected += rejected_steps(result)
        error = max(error, norm(solution(result) - target, Inf) / max(1, norm(target, Inf)))
        evaluate!(u, H, HC.DoubleDouble.ComplexDF64.(solution(result)), 0.0)
        residual = max(residual, Float64(norm(u, Inf)))
    end
    @printf(
        "%s\t%s\t%.3f\t%.3f\t%d/%d\t%.3e\t%.3e\t%d\t%d\n",
        case.name,
        string(compile),
        minimum(times) * 1e6,
        sort(times)[cld(samples, 2)] * 1e6,
        successes,
        length(case.starts),
        error,
        residual,
        accepted,
        rejected
    )
    flush(stdout)
    successes == length(case.starts) && error <= 1e-7
end

function main()
    println("# source: ", pathof(HC))
    println(
        "case\tcompile\tmin_us/path\tmedian_us/path\tsuccesses\tendpoint_error\tresidual\taccepted\trejected",
    )
    cases = tracking_cases()
    valid = true
    for compile in (false, :mixed), case in cases
        valid &= benchmark_case(case, compile)
    end
    valid || error(
        "A benchmark path failed or missed its expected endpoint; timings are not comparable.",
    )
end
main()
