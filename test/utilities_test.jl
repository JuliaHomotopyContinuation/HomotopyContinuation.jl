@testset "Utilities" begin
    @testset "UniquePoints" begin
        Random.seed!(1234)
        X = [randn(ComplexF64, 10) for _ = 1:2_000]
        indices = unique!(rand(1:2000, 20))

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

        # Test many points with nearly indentical distance to the inserted point
        points = shuffle!([[cis(k/100*2π)] for k=0:99])
        data = Utilities.UniquePoints(points)
        @test Utilities.iscontained(data, [0.0im]) == false
    end

    @testset "Polynomials" begin
        @polyvar p a b c x y z
        e = [p + 1]
        f = [a * b * c ]
        g = [x+y, y + z, x + z]
        @test Utilities.expand(f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2]
        @test Utilities.expand(e ∘ f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2 + 1]

        @test Utilities.validate(f ∘ g) == true
        @test Utilities.validate(f ∘ g, parameters=[z]) == true
        @test Utilities.validate(f ∘ g, parameters=[c]) == false
        @test Utilities.validate(g ∘ f) == false

        @test Utilities.nvariables(f ∘ g) == 3
        @test Utilities.nvariables(f ∘ g, parameters=[z]) == 2

        @polyvar x y z
        @test Utilities.ishomogenous([x^2+y^2+x*y, x^5])
        @test Utilities.ishomogenous([x^2+y^2+x*y, x^4+1]) == false

        #test weighted degree
        @test Utilities.ishomogenous(x^3+x*y, [(x, 2), (y, 4)])
        @test Utilities.homogenize(x+x*y, [(x, 2), (y, 4)], z) == x*z^4+x*y

        @test Utilities.ishomogenous(Utilities.homogenize([x^2+y^2+x*y, x^4+1]))
        @test Utilities.ishomogenous([x^2+z^2 + y, x^4+z^4], [x,z]) == false
        @test Utilities.ishomogenous(Utilities.homogenize([x^2+z^2*y, x^4+z^4*y], [x,z]), [x,z]) == true

        @test Utilities.ishomogenous(f ∘ g) == true
        h = [a * b * c  + p * c^3] ∘ [x+y, y + z, x + z]
        @test Utilities.ishomogenous(h, parameters=[p]) == true
        h2 = [a * b * c  + p] ∘ [x+y, y + z, x + z]
        @test Utilities.ishomogenous(h2, parameters=[p]) == false
    end

    @testset "Misc" begin
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
end
