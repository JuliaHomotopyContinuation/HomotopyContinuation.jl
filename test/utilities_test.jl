const Utils = HomotopyContinuation.Utilities

@testset "Utilities" begin
    PolyImpl.@polyvar x y z

    @test Utils.totaldegree([x^2+3, y^2*x^2-2]) == [x^2-1, y^4-1]
    @test Utils.totaldegree([x^2+3z^2, y^2*x^2-2z^4]) == [y^2-x^2, z^4-x^4]

    @test Utils.ishomogenous(x^2+y^2+x*y)
    @test Utils.ishomogenous([x^2+y^2+x*y, x^5])
    @test Utils.ishomogenous([x^2+y^2+x*y, x^4+1]) == false

    @test Utils.ishomogenous(Utils.homogenize([x^2+y^2+x*y, x^4+1]))
    @test Utils.nvariables(Utils.homogenize([x^2+y^2+x*y, x^4+1])) == 3
end
