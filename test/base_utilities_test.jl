using HomotopyContinuation
using Base.Test

@testset "utilities" begin

    @testset "affine/projective" begin
        @test projective([1, 2, 3]) == [1, 1, 2, 3]
        @test affine(projective([1, 2, 3])) == [1, 2, 3]
    end

end