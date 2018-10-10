@testset "Utilities" begin
    @testset "Div" begin
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

    @testset "UniquePoints" begin
        Random.seed!(1234)
        X = [randn(ComplexF64, 10) for _ = 1:2_000]
        indices = unique!(rand(1:2000, 100))

        data = Utilities.UniquePoints(X)
        @test length(data) == 2_000

        for i ∈ indices
            @test Utilities.iscontained(data, X[i])
            @test Utilities.iscontained(data, X[i], Val(true)) == i
            @test Utilities.iscontained(data, X[i] .+ 1e-4) == false
            @test Utilities.iscontained(data, X[i] .+ 1e-9, Val(true)) == i
            @test Utilities.iscontained(data, X[i] .+ 1e-9) == true
            @test Utilities.add!(data, X[i]) == false
            @test Utilities.add!(data, X[i], Val(true)) == i
            @test data[i] == X[i]
            @test Utilities.add!(data, X[i] .+ 1e-4) == true
            @test Utilities.add!(data, X[i] .- 1e-4, Val(true)) == -1
        end
    end
end
