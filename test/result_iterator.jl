@testset "ResultIterator" begin
    @testset "Two Plane Curves" begin
        @var x y
        # define the polynomials
        f₁ = y-x^2
        f₂ = y-x^3
        F = [f₁, f₂]
        tsi_polyhedral = solve(F; iterator_only=true, start_system = :polyhedral)
        tsi_total_degree = solve(F; iterator_only=true, start_system = :total_degree)

        @test nsolutions(tsi_polyhedral) == 3
        @test nsolutions(tsi_total_degree) == 3
        @test length(tsi_polyhedral) == 3
        @test length(tsi_total_degree) == 6
        @test bitmask(tsi_total_degree,isfinite) == Vector{Bool}([1,1,1,0,0,0]) #! Why is this not deterministic?!
        
    end

    @testset "Two Plane Curves - with bitmask filtering" begin
        @var x y
        # define the polynomials
        f₁ = y-x^2
        f₂ = y-x^3
        F = [f₁, f₂]
        tsi_polyhedral = solve(F; iterator_only=true, start_system = :polyhedral)
        tsi_total_degree = solve(F; iterator_only=true, start_system = :total_degree)
        #tsi_non_singular = HomotopyContinuation.bitmask_filter(tsi_total_degree,x->!is_singular(x))
        #tsi_filtered = HomotopyContinuation.bitmask_filter(tsi_total_degree)
    end
end