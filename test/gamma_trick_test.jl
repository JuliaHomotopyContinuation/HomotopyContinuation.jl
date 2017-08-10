@testset "gamma trick" begin
    
    @TP.polyvar x y a
    f = [x^2+(3.0+0.0im)*y]
    g = [y^2-(1.0+0.0im)*x]

    H = GammaTrickHomotopy(g, f, 1.0+0im)
    @test typeof(H)<:GammaTrickHomotopy{Complex128}

    @test evaluate(H, [2.0+0im, 1.0], 1.0) == [-1.0+0im]
    @test evaluate(H, [2.0+0im, 1.0], 0.0) == [7.0+0im]

    J_H = jacobian(H)

    @test J_H([1.0+0im, 1.0], 0.0) == [2 3+0im]
    @test J_H([1.0+0im, 1.0], 1.0) == [-1 2+0im]

    ∂H∂t = dt(H)

    ## time derivate is independent of t
    @test ∂H∂t([1.0+0im,1.0], 0.23) == [-4+0im]

    @test degrees(H) == [2]
    @test startsystem(H) == g
    @test targetsystem(H) == f


    @testset "homogenize" begin
        @TP.polyvar x y z
        f = [x^2+(3.0+0.0im)*y, 2y-(3.0+0.0im)*x]
        g = [y^2-(1.0+0.0im)*x, (3.0+0.0im)*x + y]

        H = GammaTrickHomotopy(g, f)

        K = homogenize(H, z)
        
        @test nvars(K) == 3
        @test H.γ == K.γ
    end

    @testset "constructor" begin
        @TP.polyvar x y a
        f = [x^2+(3.0+0.0im)*y]
        g = [y^2-(1.0+0.0im)*x]

        H = GammaTrickHomotopy(g, f)
        @test norm(H.γ) ≈ 1.0 atol=1e-8

        H = GammaTrickHomotopy(g, f, 212)
        K = GammaTrickHomotopy(g, f, 212)
        @test H.γ ≈ K.γ

        H = GammaTrickHomotopy(g, f, 1.3+4im)
        @test H.γ ≈ 1.3+4im

        @test_throws ErrorException GammaTrickHomotopy(g, [MP.polynomial((1.0+0.0im)a)], 1.0+0im)
        @test_throws ErrorException GammaTrickHomotopy([x^2+(3.0+0.0im)y, y^2-(1.0+0.0im)x], [x^2+(3.0+0.0im)y], 1.0+0im)
    end

end