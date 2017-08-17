@testset "GammaTrickHomotopy" begin

    PolyImpl.@polyvar x y

    H = GammaTrickHomotopy([y^2-x], [x^2+3y], 1.0+0im)
    @test typeof(H)<:GammaTrickHomotopy{Int64}

    @test evaluate(H, [2.0+0im, 1.0], 1.0) == [-1.0+0im]
    @test H([2.0+0im, 1.0], 1.0) == [-1.0+0im]
    @test evaluate(H, [2.0+0im, 1.0], 0.0) == [7.0+0im]

    J_H = differentiate(H)

    @test J_H([1.0+0im, 1.0], 0.0) == [2 3+0im]
    @test J_H([1.0+0im, 1.0], 1.0) == [-1 2+0im]

    ∂H∂t = dt(H)

    ## time derivate is independent of t
    @test ∂H∂t([1.0+0im,1.0], 0.23) == [-4+0im]

    @test degrees(H) == [2]
    @test startsystem(H) == PolySystem([y^2-x])
    @test targetsystem(H) == PolySystem([x^2+3.0y])


    f = [x^2+3y, 2y-3x]
    g = [y^2-x, 3x + y]
    H = GammaTrickHomotopy(g, f)
    @test homogenized(H) == false
    K = homogenize(H)
    @test ishomogenous(H) == false
    @test homogenized(K)
    @test ishomogenous(K)
    @test nvariables(K) == 3
    @test H.γ == K.γ

    H = GammaTrickHomotopy([y^2-x], [x^2+3y])
    @test norm(H.γ) ≈ 1.0 atol=1e-8

    H = GammaTrickHomotopy(g, f, 212)
    K = GammaTrickHomotopy(g, f, 212)
    @test H.γ ≈ K.γ

    H = GammaTrickHomotopy(g, f, 1.3+4im)
    @test H.γ ≈ 1.3+4im

    PolyImpl.@polyvar a

    @test_throws ErrorException GammaTrickHomotopy([x^2+3y], [a], 1.0+0im)
    @test_throws ErrorException GammaTrickHomotopy([x^2+3y, y^2-x], [x^2+3y], 1.0+0im)
end
