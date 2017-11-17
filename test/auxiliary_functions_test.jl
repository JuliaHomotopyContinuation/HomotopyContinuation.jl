@testset "isolated zeros" begin
    @polyvar x y
    f = [x^2-1, y^2-1]
    H,s = totaldegree(GeodesicOnTheSphere, f)
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(real_solutions(Sol)) == 4
    @test length(singular_solutions(Sol)) == 0
    @test length(isolated_solutions(Sol)) == 4
end


@testset "singular zeros" begin
    @polyvar x
    f = (x+1)^2
    H,s = totaldegree(GeodesicOnTheSphere, [f])
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(real_solutions(Sol)) == 2
    @test length(singular_solutions(Sol)) == 2
    @test length(isolated_solutions(Sol)) == 0
end


@testset "solutions at at infinity" begin
    @polyvar x y
    f = [(x-1)*(y+2), (x-2)*(y-2)]
    H,s = totaldegree(StraightLineHomotopy, f)
    cfg=PolynomialHomotopyConfig(H)
    Sol=solve(H,s);

    @test length(real_solutions(Sol)) == 4
    @test length(isolated_solutions(Sol)) == 2
    @test length(solutions_at_infinity(Sol)) == 2
end
