@testset "Linear Spaces" begin
    @testset "AffineSubspaceHomotopy" begin
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)

        A = rand_affine_subspace(4; codim = 2)

        # Commpute witness set
        L = A(x, Extrinsic)
        H, S = total_degree_homotopy(F ∩ L)
        res = track.(Tracker(H), S)
        W = solution.(res)

        B = rand_affine_subspace(4, codim = 2)

        graff_tracker = Tracker(AffineSubspaceHomotopy(F, A, B))
        graff_result = track.(graff_tracker, W)
        @test all(is_success, graff_result)
    end
end
