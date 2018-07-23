@testset "Utilities" begin
    @polyvar x y z
    @test Utilities.ishomogenous(x^2+y^2+x*y)
    @test Utilities.ishomogenous([x^2+y^2+x*y, x^5])
    @test Utilities.ishomogenous([x^2+y^2+x*y, x^4+1]) == false

    @test Utilities.ishomogenous(Utilities.homogenize([x^2+y^2+x*y, x^4+1]))
    @test Utilities.nvariables(Utilities.homogenize([x^2+y^2+x*y, x^4+1])) == 3

    @test Utilities.ishomogenous([x^2+z^2 + y, x^4+z^4], [x,z]) == false
    @test Utilities.ishomogenous(Utilities.homogenize([x^2+z^2*y, x^4+z^4*y], [x,z]), [x,z]) == true

    A = rand(Complex{Float64}, 12, 12)
    b = rand(Complex{Float64}, 12)
    C, d = copy(A), copy(b)
    @test norm(Utilities.solve!(C, d) - A \ b) < 1e-10

    A = rand(Complex{Float64}, 15, 12)
    b = rand(Complex{Float64}, 15)
    C, d = copy(A), copy(b)
    Utilities.solve!(C, d)
    @test norm(d[1:12] - A \ b) < 1e-10
end
