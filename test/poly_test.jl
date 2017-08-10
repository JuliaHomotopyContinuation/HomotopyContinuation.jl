@testset "Poly" begin
    @TP.polyvar x y z

    f = 3x^2-2y^2*x-y+(3.0+0.0im)
    @test evaluate(f, [2, 1]) ≈ 3*4-2*1*2-1+3.0+0.0im
    
    @test gradient(f) == [ 6x - 2y^2, -1 - 4x*y]

    @test is_homogenous(f) == false

    # currently there's a bug in TP so that the equality fails (the result is correct)
    @test_broken homogenize(f,z) == (3.0+0.0im)*z^3 - y*z^2 + 3x^2*z - 2x*y^2

    @test is_homogenous(homogenize(f,z)) == true

    @test deg(f) == 3

    @test nvars(f) == 2

    f, g = promote(3.0x^2 + 2.0x * y - 1.0 * y^2, (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^2)    
    @test weyl_dot(f, g) == 3 * conj(-2.5 + 2im) - 2 * 3 / 2 - 4
    @test weyl_dot(f, f) ==  9.0 + 4.0  / 2 + 1.0 
end