@testset "TotalDegree" begin
    @test HC.bezout_number([2 1; 2 1; 1 1], [2, 1]) == 8
    @test HC.bezout_number([2 1 0; 1 1 1; 0 1 1], [1, 1, 1]) == 5
end
