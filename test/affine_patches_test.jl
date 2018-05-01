@testset "AffinePatches" begin
    patch = AffinePatches.OrthogonalPatch()
    @test patch isa AffinePatches.AbstractLocalAffinePatch
    @test patch isa AffinePatches.AbstractAffinePatch

    x = rand(Complex{Float64}, 10)
    v = AffinePatches.init_patch(patch, x)
    @test v == x
    AffinePatches.precondition!(v, x, patch)
    @test norm(x) ≈ 1.0
    @test v == x

    y = rand(Complex{Float64}, 10)
    AffinePatches.update_patch!(v, patch, y)
    @test v == y
end
