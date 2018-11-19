@testset "PathTracking" begin
    @testset "General" begin
        F = equations(katsura(5))

        P, start_sols = Problems.problem_startsolutions(Input.TotalDegree(F))
        x₁ = Problems.embed(P, first(start_sols))
        H = Homotopies.PatchedHomotopy(P.homotopy, AffinePatches.OrthogonalPatch(), x₁)
        # test construction
        t1 = PathTracking.PathTracker(H, x₁, 1.0, 0.1)

        @test_nowarn PathTracking.currΔt(t1)

        rtol = PathTracking.refinement_tol(t1)
        @test_nowarn PathTracking.set_refinement_tol!(t1, 5e-5)
        @test PathTracking.refinement_tol(t1) == 5e-5
        PathTracking.set_refinement_tol!(t1, rtol)

        rmaxiter = PathTracking.refinement_maxiters(t1)
        @test_nowarn PathTracking.set_refinement_maxiters!(t1, 11)
        @test PathTracking.refinement_maxiters(t1) == 11
        PathTracking.set_refinement_maxiters!(t1, rmaxiter)

        @test t1 isa PathTracking.PathTracker
        @test length(t1.x) == 7

        t2 = PathTracking.PathTracker(H, x₁, 1.0, 0.1)
        @test t2 isa PathTracking.PathTracker
        @test length(t2.x) == 7

        t3 = PathTracking.PathTracker(H, x₁, 1.0, 0.1, predictor=Predictors.Euler())
        @test t3.predictor_corrector.predictor isa Predictors.Euler


        PathTracking.setup!(t1, Problems.embed(P, first(start_sols)), 1.0, 0.4)
        @test PathTracking.currstatus(t1) == :ok
        @test PathTracking.currt(t1) == 1.0

        PathTracking.setup!(t1, Problems.embed(P, first(start_sols)), 0.5, 0.4)
        @test PathTracking.currstatus(t1) == :invalid_startvalue
        @test PathTracking.currt(t1) == 0.5

        R = PathTracking.track(t1, Problems.embed(P, first(start_sols)), 1.0, 0.0)
        @test R isa PathTracking.PathTrackerResult
        @test R.returncode == :success
        @test R.res < 1e-7
        @test_nowarn show(devnull, R)
        @test_nowarn TreeViews.treelabel(devnull, R, MIME"application/prs.juno.inline"())

        out = Problems.embed(P, first(start_sols))
        retcode = PathTracking.track!(out, t1, Problems.embed(P, first(start_sols)), 1.0, 0.0)
        @test retcode == :success
        @test out == R.x
    end

    @testset "LinearSystem" begin
        A = rand(3, 3)
        b = rand(3)

        @polyvar x y z
        F = A * [x, y, z] - b

        P, sols = Problems.problem_startsolutions(Input.TotalDegree(F))

        start_sols = map(s -> Problems.embed(P,s), sols)

        s = start_sols[1]
        H = Homotopies.PatchedHomotopy(P.homotopy, AffinePatches.RandomPatch(), s)
        # intilizae tracker with random start and target to see whether it is properly resetetted
        tracker = PathTracking.PathTracker(H, s, rand(), rand())

        result = PathTracking.track(tracker, s, 1.0, 0.0)
        @test result.returncode == :success
        @test result.t == 0.0
        x = result.x
        @test norm(affine(x) - A \ b) < 1e-6

        x_inter = copy(s)
        retcode = PathTracking.track!(x_inter, tracker, s, 1.0, 0.1)
        @test retcode == :success
        x_final = zero(x_inter)
        retcode = PathTracking.track!(x_final, tracker, x_inter, 0.1, 0.0)
        @test retcode == :success
        tracker
        @test PathTracking.curriters(tracker) < 3
        x = PathTracking.currx(tracker)
        @test norm(affine(x) - A \ b) < 1e-6
    end

    @testset "fixedpatch" begin
        F = equations(katsura(5))
        P, start_sols = Problems.problem_startsolutions(Input.TotalDegree(F))
        x1 = Problems.embed(P, first(start_sols))
        patch = AffinePatches.state(AffinePatches.OrthogonalPatch(), x1)
        fixedpatch = AffinePatches.state(AffinePatches.FixedPatch(), x1)
        H = Homotopies.PatchedHomotopy(P.homotopy, patch)
        FH = Homotopies.PatchedHomotopy(P.homotopy, fixedpatch)
        # test construction
        tracker = PathTracking.PathTracker(H, x1, 1.0, 0.1)
        fixedtracker = PathTracking.PathTracker(FH, x1, 1.0, 0.1)

        r1 = PathTracking.track(fixedtracker, x1, 1.0, 0.1)
        v1 = copy(fixedpatch.v_conj)
        PathTracking.track!(fixedtracker, r1.x, 0.1, 0.05, precondition=false)
        @test v1 == fixedpatch.v_conj
        PathTracking.track!(fixedtracker, PathTracking.currx(tracker), 0.05, 0.01, precondition=false)
        @test v1 == fixedpatch.v_conj
        PathTracking.track!(fixedtracker, r1.x, 0.1, 0.05, precondition=true)
        @test v1 != fixedpatch.v_conj

        r1 = PathTracking.track(tracker, x1, 1.0, 0.1)
        v1 = copy(patch.v_conj)
        PathTracking.track!(tracker, r1.x, 0.1, 0.0)
        @test v1 != patch.v_conj
    end
end
