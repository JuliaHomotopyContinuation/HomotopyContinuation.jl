@testset "Linear Spaces" begin
    @testset "LinearSubspace" begin
        A = LinearSubspace([1 0 3; 2 1 3], [5, -2])
        @test dim(A) == dim(intrinsic(A)) == dim(extrinsic(A)) == 1
        @test codim(A) == codim(intrinsic(A)) == codim(extrinsic(A)) == 2
        @test ambient_dim(A) == 3
        @test startswith(sprint(show, A), "1-dim. affine linear subspace")
        @test identity.(A) == A

        @test intrinsic(A) isa IntrinsicDescription
        @test startswith(sprint(show, intrinsic(A)), "IntrinsicDescription")
        @test extrinsic(A) isa ExtrinsicDescription
        @test startswith(sprint(show, extrinsic(A)), "ExtrinsicDescription")

        @test identity.(extrinsic(A)) == extrinsic(A)
        @test identity.(intrinsic(A)) == intrinsic(A)
        @test identity.(Intrinsic) == Intrinsic

        u = rand(1)
        x = A(u, Intrinsic)
        @test coord_change(A, Extrinsic, Intrinsic, x) ≈ u rtol = 1e-12
        @test coord_change(A, Intrinsic, Extrinsic, u) ≈ x rtol = 1e-12
        @test norm(A(x, Extrinsic)) ≈ 0 atol = 1e-14
        @test norm(ExtrinsicDescription(intrinsic(A))(x)) ≈ 0.0 atol = 1e-12

        B = rand_subspace(3; codim = 2)
        @test B != A
        copy!(B, A)
        @test extrinsic(B) == extrinsic(A)
        @test intrinsic(B) == intrinsic(A)
        @test B == A
        @test copy(A) == A
        @test copy(A) !== A

        C = rand_subspace(3, dim = 1, real = true)
        @test geodesic_distance(A, C) > 0
        γ = geodesic(A, C)
        @test size(γ(1)) == size(intrinsic(A).Y)

        A2 = translate(A, [1, 1], Extrinsic)
        A3 = LinearSubspace(extrinsic(A).A, extrinsic(A).b + [1, 1])
        @test A2.intrinsic.Y ≈ A3.intrinsic.Y

        # rand subspace through a point
        x = randn(ComplexF64, 5)
        L = rand_subspace(x; dim = 2)
        @test norm(L(x)) ≈ 0 atol = 1e-8
        L2 = rand_subspace(x; dim = 2, affine = false)
        @test is_linear(L2)
        @test norm(L2(x)) ≈ 0 atol = 1e-8
    end

    @testset "LinearSubspaceHomotopy" begin
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)

        A = rand_subspace(4; codim = 2)
        # Commpute witness set
        L = A(x, Extrinsic)
        res = track.(total_degree(F ∩ L)...)
        W = solution.(res)

        B = rand_subspace(4, codim = 2)

        graff_tracker = Tracker(
            LinearSubspaceHomotopy(F, A, B),
            options = (automatic_differentiation = 3,),
        )
        graff_result = track.(graff_tracker, W)
        @test all(is_success, graff_result)

        C = rand_subspace(4, codim = 2)
        copy_A = copy(A)
        set_subspaces!(graff_tracker.homotopy, B, C)
        @test A == copy_A
        @test all(is_success, track.(graff_tracker, solution.(graff_result)))

        graff_path_tracker = EndgameTracker(LinearSubspaceHomotopy(F, A, B))
        graff_path_result = track.(graff_path_tracker, W)
        @test all(is_success, graff_path_result)

        @test all(isapprox.(
            solution.(graff_path_result),
            solution.(graff_result),
            rtol = 1e-12,
        ))
    end
end
