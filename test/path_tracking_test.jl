@testset "PathTracking" begin
    @testset "General" begin
        F = equations(katsura(5))
        # test construction
        t1, start_sols = pathtracker_startsolutions(F, patch=OrthogonalPatch())

        @test_nowarn currΔt(t1)

        raccuracy = refinement_accuracy(t1)
        @test_nowarn set_refinement_accuracy!(t1, 5e-5)
        @test refinement_accuracy(t1) == 5e-5
        set_refinement_accuracy!(t1, raccuracy)

        rmaxiter = max_refinement_iters(t1)
        @test_nowarn set_max_refinement_iters!(t1, 11)
        @test max_refinement_iters(t1) == 11
        set_max_refinement_iters!(t1, rmaxiter)

        @test t1 isa PathTracker
        @test length(currx(t1)) == 7

        t3, start_sols = pathtracker_startsolutions(F, patch=OrthogonalPatch(), predictor=Euler())
        @test t3.predictor isa Euler

        setup!(t1, first(start_sols), 1.0, 0.4)
        @test currstatus(t1) == PathTrackerStatus.tracking
        @test currt(t1) == 1.0

        setup!(t1, first(start_sols), 0.5, 0.4)
        @test currstatus(t1) == PathTrackerStatus.terminated_invalid_startvalue
        @test currt(t1) == 0.5

        R = track(t1, first(start_sols), 1.0, 0.0)
        @test R isa PathTrackerResult
        @test R.returncode == PathTrackerStatus.success
        @test R.accuracy < 1e-7
        @test_nowarn show(devnull, R)

        out = [1; first(start_sols)]
        retcode = track!(out, t1, first(start_sols), 1.0, 0.0)
        @test retcode == PathTrackerStatus.success
        @test out == R.x
    end

    @testset "Affine tracking" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im], parameters=p, p₁=[1, 0], p₀=[2, 4],
                    affine=true)
        @test affine_tracking(tracker) == true
        res = track(tracker, starts[1], 1.0, 0.0)
        @test res.returncode == PathTrackerStatus.success
        @test isa(res.x, Vector{ComplexF64})
        @test length(res.x) == 2

        x = @SVector [1.0, 1.0 + 0.0*im]
        tracker, starts = pathtracker_startsolutions(F, x, parameters=p, p₁=[1, 0], p₀=[2, 4], affine=true)
        @test length(starts) == 1
        res = track(tracker, starts[1], 1.0, 0.0)
        @test isa(res.x, Vector{ComplexF64})
        @test length(res.x) == 2
    end

    @testset "Allocations" begin
        F = equations(katsura(5))
        tracker, start_sols = pathtracker_startsolutions(F)
        s = first(start_sols)
        @allocated track!(tracker, s, 1.0, 0.1)
        @test (@allocated track!(tracker, s, 1.0, 0.1)) == 0
    end

    @testset "LinearSystem" begin
        A = rand(3, 3)
        b = rand(3)

        @polyvar x[1:3]
        F = A * x - b

        tracker, start_sols = pathtracker_startsolutions(F, patch=RandomPatch(), max_steps=10)
        s = first(start_sols)
        result = track(tracker, s, 1.0, 0.0)
        @test result.returncode == PathTrackerStatus.success
        @test result.t == 0.0
        x = result.x
        @test norm(ProjectiveVectors.affine_chart(x) - A \ b) < 1e-6

        x_inter = [1;copy(s)]
        retcode = track!(x_inter, tracker, s, 1.0, 0.1)
        @test retcode == PathTrackerStatus.success
        x_final = zero(x_inter)
        retcode = track!(x_final, tracker, x_inter, 0.1, 0.0)
        @test retcode == PathTrackerStatus.success
        @test curriters(tracker) < 3
        x = currx(tracker)
        @test norm(ProjectiveVectors.affine_chart(x) - A \ b) < 1e-6
    end

    @testset "FixedPatch continuation" begin
        F = equations(katsura(5))
        # test construction
        tracker, start_sols = pathtracker_startsolutions(F, patch=OrthogonalPatch())
        patch = tracker.state.patch
        fixedtracker = pathtracker(F, patch=FixedPatch())
        fixedpatch = fixedtracker.state.patch
        x1 = first(start_sols)

        r1 = track(fixedtracker, x1, 1.0, 0.1)
        v1 = copy(fixedpatch.v̄)
        track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch=false)
        @test v1 == fixedpatch.v̄
        track!(fixedtracker, currx(tracker), 0.05, 0.01, setup_patch=false)
        @test v1 == fixedpatch.v̄
        track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch=true)
        @test v1 != fixedpatch.v̄

        r1 = track(tracker, x1, 1.0, 0.1)
        v1 = copy(patch.v̄)
        track!(tracker, r1.x, 0.1, 0.0)
        @test v1 != patch.v̄
    end

    @testset "pathtracker_startsolutions" begin
        @polyvar  x y

        tracker, S = pathtracker_startsolutions([x^2-y^2+4, x + y - 3])
        result = track(tracker, first(S), 1.0, 0.0)
        @test result.returncode == PathTrackerStatus.success

        # test type promotion of startsolutions
        tracker, S = pathtracker_startsolutions([x^2-y^2+4, x + y - 3])
        result = track(tracker, [1, 1], 1.0, 0.0)
        @test result.returncode == PathTrackerStatus.success
    end

    @testset "iterator" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]
        s = [1.0, 1.0 + 0.0*im]
        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4])
        typeof(first(iterator(tracker, s, 1.0, 0.0))) ==
            Tuple{Vector{ComplexF64}, Float64}
        typeof(first(iterator(tracker, s, 1.0, 0.0; affine=false))) ==
            Tuple{ProjectiveVectors.PVector{ComplexF64, 1}, Float64}

        @test max_step_size(tracker) == Inf
        set_max_step_size!(tracker, 0.01)
        @test max_step_size(tracker) == 0.01

        length(collect(iterator(tracker, s, 1.0, 0.0))) == 101


        # Test iteration protocol
        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4])
        setup!(tracker, s, 1.0, 0.0)
        for tr in tracker
            # nothing
        end
        @test tracker.state.status == PathTrackerStatus.success
    end

    @testset "Deprecation warnings" begin
        F = equations(katsura(5))
        # test deprecation mechanism
        tracker = pathtracker(F, tol=1e-8)
        @test accuracy(tracker) == 1e-8
    end
end
