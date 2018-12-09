@testset "PathTracking" begin
    @testset "General" begin
        F = equations(katsura(5))
        # test construction
        t1, start_sols = PathTracking.pathtracker_startsolutions(F, patch=AffinePatches.OrthogonalPatch())

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
        @test length(PathTracking.currx(t1)) == 7

        t3, start_sols = PathTracking.pathtracker_startsolutions(F, patch=AffinePatches.OrthogonalPatch(), predictor=Predictors.Euler())
        @test t3.predictor isa Predictors.Euler

        PathTracking.setup!(t1, first(start_sols), 1.0, 0.4)
        @test PathTracking.currstatus(t1) == PathTracking.Status.tracking
        @test PathTracking.currt(t1) == 1.0

        PathTracking.setup!(t1, first(start_sols), 0.5, 0.4)
        @test PathTracking.currstatus(t1) == PathTracking.Status.terminated_invalid_startvalue
        @test PathTracking.currt(t1) == 0.5

        R = PathTracking.track(t1, first(start_sols), 1.0, 0.0)
        @test R isa PathTracking.PathTrackerResult
        @test R.returncode == PathTracking.Status.success
        @test R.accuracy < 1e-7
        @test_nowarn show(devnull, R)

        out = [1; first(start_sols)]
        retcode = PathTracking.track!(out, t1, first(start_sols), 1.0, 0.0)
        @test retcode == PathTracking.Status.success
        @test out == R.x
    end

    @testset "LinearSystem" begin
        A = rand(3, 3)
        b = rand(3)

        @polyvar x[1:3]
        F = A * x - b

        tracker, start_sols = pathtracker_startsolutions(F, patch=AffinePatches.RandomPatch())
        s = first(start_sols)
        result = PathTracking.track(tracker, s, 1.0, 0.0)
        @test result.returncode == PathTracking.Status.success
        @test result.t == 0.0
        x = result.x
        @test norm(ProjectiveVectors.affine(x) - A \ b) < 1e-6

        x_inter = [1;copy(s)]
        retcode = PathTracking.track!(x_inter, tracker, s, 1.0, 0.1)
        @test retcode == PathTracking.Status.success
        x_final = zero(x_inter)
        retcode = PathTracking.track!(x_final, tracker, x_inter, 0.1, 0.0)
        @test retcode == PathTracking.Status.success
        tracker
        @test PathTracking.curriters(tracker) < 3
        x = PathTracking.currx(tracker)
        @test norm(ProjectiveVectors.affine(x) - A \ b) < 1e-6
    end

    @testset "FixedPatch continuation" begin
        F = equations(katsura(5))
        # test construction
        tracker, start_sols = PathTracking.pathtracker_startsolutions(F, patch=AffinePatches.OrthogonalPatch())
        patch = tracker.state.patch
        fixedtracker = PathTracking.pathtracker(F, patch=AffinePatches.FixedPatch())
        fixedpatch = fixedtracker.state.patch
        x1 = first(start_sols)

        r1 = PathTracking.track(fixedtracker, x1, 1.0, 0.1)
        v1 = copy(fixedpatch.v̄)
        PathTracking.track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch=false)
        @test v1 == fixedpatch.v̄
        PathTracking.track!(fixedtracker, PathTracking.currx(tracker), 0.05, 0.01, setup_patch=false)
        @test v1 == fixedpatch.v̄
        PathTracking.track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch=true)
        @test v1 != fixedpatch.v̄

        r1 = PathTracking.track(tracker, x1, 1.0, 0.1)
        v1 = copy(patch.v̄)
        PathTracking.track!(tracker, r1.x, 0.1, 0.0)
        @test v1 != patch.v̄
    end
end
