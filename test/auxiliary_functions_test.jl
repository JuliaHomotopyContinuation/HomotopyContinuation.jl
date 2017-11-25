@testset "isolated zeros" begin
    PolyImpl.@polyvar x y
    f = [x^2-1, y^2-1]
    H,s = totaldegree(GeodesicOnTheSphere, f)
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(solutions(Sol, only_real = true)) == 4
    @test length(solutions(Sol, singular = false)) == 4
end


@testset "singular zeros" begin
    PolyImpl.@polyvar x
    f = (x+1)^2
    H,s = totaldegree(GeodesicOnTheSphere, [f])
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(solutions(Sol, only_real = true)) == 2
    @test length(solutions(Sol, singular = true)) == 2
end


@testset "solutions at at infinity" begin
    PolyImpl.@polyvar x y
    f = [(x-1)*(y+2), (x-2)*(y-2)]
    H,s = totaldegree(StraightLineHomotopy, f)
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(solutions(Sol, at_infinity = false)) == 2
end
