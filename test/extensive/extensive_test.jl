@testset "Lines on a Quintic surface in 3-space" begin
    @testset "Lines on a quintic surface in 3-space" begin
        sys, q₀ = fano_quintic()
        @time res =
            solve(sys; gamma = gamma, target_parameters = q₀, start_system = :total_degree)
        @test nsolutions(res) == 2875

        @time poly_res = solve(sys; target_parameters = q₀, start_system = :polyhedral)
        @test nsolutions(poly_res) = 2875
    end

end
