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
        γ1 = γ(1)
        @test ambient_dim(γ1) == ambient_dim(A)
        @test dim(γ1) == dim(A)           

        A2 = translate(A, [1, 1], Extrinsic)
        A3 = LinearSubspace(extrinsic(A).A, extrinsic(A).b + [1, 1])
        @test A2.intrinsic.X ≈ A3.intrinsic.X

        # rand subspace through a point
        x = randn(ComplexF64, 5)
        L = rand_subspace(x; dim = 2)
        @test norm(L(x)) ≈ 0 atol = 1e-8
        L2 = rand_subspace(x; dim = 2, affine = false)
        @test is_linear(L2)
        @test norm(L2(x)) ≈ 0 atol = 1e-8
    end

    @testset "Intersect subspaces" begin
        L₁ = rand_subspace(7; codim = 2)
        L₂ = rand_subspace(7; codim = 3)
        L₃ = L₁ ∩ L₂
        @test codim(L₃) == 5
        @test dim(L₃) == 2
        @test ambient_dim(L₃) == 7
        E₃ = extrinsic(L₃)
        @test norm(L₁(E₃.A \ E₃.b, Extrinsic)) ≈ 0 atol = 1e-12
        @test norm(L₂(E₃.A \ E₃.b, Extrinsic)) ≈ 0 atol = 1e-12
        @test norm(L₃(E₃.A \ E₃.b, Extrinsic)) ≈ 0 atol = 1e-12
    end

    @testset "IntrinsicSubspaceHomotopy" begin
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
            IntrinsicSubspaceHomotopy(F, A, B),
            options = (automatic_differentiation = 3,),
        )
        graff_result = track.(graff_tracker, W)
        @test all(is_success, graff_result)

        C = rand_subspace(4, codim = 2)
        copy_A = copy(A)
        set_subspaces!(graff_tracker.homotopy, B, C)
        @test A == copy_A
        @test all(is_success, track.(graff_tracker, solution.(graff_result)))

        target_parameters!(graff_tracker.homotopy, C)
        @test A == copy_A
        @test all(is_success, track.(graff_tracker, solution.(graff_result)))

        graff_path_tracker = EndgameTracker(IntrinsicSubspaceHomotopy(F, A, B))
        graff_path_result = track.(graff_path_tracker, W)
        @test all(is_success, graff_path_result)

        @test all(
            isapprox.(solution.(graff_path_result), solution.(graff_result), rtol = 1e-12),
        )
    end

    @testset "SubspaceHomotopy between perpendicular spaces" begin
        @var x, y, z
        p = (x * y - x^2) + 1 - z
        q = x^4 + x^2 - y - 1
        f = [
            p * q * (x - 3) * (x - 5)
            p * q * (y - 3) * (y - 5)
            p * (z - 3) * (z - 5)
        ]
        F = System(f, variables = [x, y, z])

        L1 = LinearSubspace(reshape([1; 0; 0.0], 1 ,3), [1.0])
        L2 = LinearSubspace(reshape([0; 1; 0.0], 1, 3), [1.0])

        W = witness_set(F, L1; show_progress = false, start_system = :total_degree)
        start = solutions(W)[1]

        H1 = ExtrinsicSubspaceHomotopy(F, L1, L2);
        H2 = IntrinsicSubspaceHomotopy(F, L1, L2);
        H3 = IntrinsicSubspaceProjectiveHomotopy(F, L1, L2);

        T1 = EndgameTracker(H1);
        T2 = EndgameTracker(H2);
        T3 = EndgameTracker(H3);

        res1 = track(T1, start)
        res2 = track(T2, start)
        res3 = track(T3, start)

        @test is_success(res1)
        @test is_success(res2)
        @test is_success(res3)
    end
end
