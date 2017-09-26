@testset "Result" begin
    a, b, c = 5, 4, 3
    PolyImpl.@polyvar sθ cθ
    f1 = cθ^2 + sθ^2 - (1.0 + 0im)
    f2 = (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2

    H, start_solutions = totaldegree(StraightLineHomotopy, [f1, f2])

    #now we can solve
    results = solve(H, start_solutions, SphericalPredictorCorrector(), report_progress=true)

    result = first(filter(issuccessfull, results))
    # check accessors
    @test norm(H(solution(result), 0.0)) ≈ 0 atol=1e-8
    @test iterations(result) > 0
    @test endgame_iterations(result) > 0
    @test algorithm(result) isa SphericalPredictorCorrector
    @test affine_solution(result) == solution(result)
    @test length(projective_solution(result)) == 3
    @test any(x -> x == startvalue(result), start_solutions)
    # currently cauchyendgame always tracks the trace and steps...
    @test length(pathtrace(result)) > 0
    @test length(first(pathtrace(result))) == 2
    @test length(pathsteps(result)) > 0
    @test length(pathsteps(result)) == length(pathtrace(result))

    # check that show doesn't throw
    @test length(string(result)) > 0

    # Now we split the it into trackpath and a cauchy endgame

    # Check result from trackpath
    HH = homogenize(H)
    start = [1; startvalue(result)]
    pathresult = trackpath(HH, start, SphericalPredictorCorrector(), 1.0, 0.1)

    @test norm(HH(HomotopyContinuation.result(pathresult), 0.1)) ≈ 0 atol=1e-8
    @test issuccessfull(pathresult)
    @test iterations(pathresult) > 0
    @test startvalue(pathresult) == start
    @test laststep(pathresult) ≈ 0.1
    @test length(pathtrace(pathresult)) == 0
    @test length(pathsteps(pathresult)) == 0

    @test length(string(pathresult)) > 0

    endgameresult = cauchyendgame(
        HH, jacobian!(HH), dt!(HH), HomotopyContinuation.result(pathresult), 0.1,
        SphericalPredictorCorrector(), geometric_series_factor=0.45)

    @test norm(HH(HomotopyContinuation.result(endgameresult), 0.0)) ≈ 0 atol=1e-8
    @test issuccessfull(endgameresult)
    @test iterations(endgameresult) > 0
    # currently cauchyendgame always tracks the trace and steps...
    @test length(pathtrace(endgameresult)) > 0
    @test length(pathsteps(endgameresult)) > 0

    @test !isnull(convergent_cluster(endgameresult))

    @test geometric_series_factor(endgameresult) == 0.45
    # test the cluster
    cc = get(convergent_cluster(endgameresult))
    @test length(clusterpoints(cc)) == 1
    @test norm(HH(cluster_convergence_point(cc), clustertime(cc))) ≈ 0 atol=1e-2
    @test length(string(cc)) > 0

    @test endgameradius(endgameresult) == clustertime(cc)
    @test length(string(endgameresult)) > 0
end
