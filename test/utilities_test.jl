@testset "Utilities" begin
    @testset "UniquePoints" begin
        Random.seed!(1234)
        X = [randn(ComplexF64, 10) for _ = 1:2_000]
        indices = HC.unique!(rand(1:2000, 20))

        data = HC.UniquePoints(X)
        @test length(data) == 2_000

        for i ∈ indices
            @test HC.iscontained(data, X[i])
            @test HC.iscontained(data, X[i], Val(true)) == i
            @test HC.iscontained(data, X[i] .+ 1e-4) == false
            @test HC.iscontained(data, X[i] .+ 1e-9, Val(true)) == i
            @test HC.iscontained(data, X[i] .+ 1e-9) == true
            @test HC.add!(data, X[i]) == false
            @test HC.add!(data, X[i], Val(true)) == i
            @test data[i] == X[i]
            @test HC.add!(data, X[i] .+ 1e-4) == true
            @test HC.add!(data, X[i] .- 1e-4, Val(true)) == -1
        end

        # Test many points with nearly indentical distance to the inserted point
        points = shuffle!([[cis(k/100*2π)] for k=0:99])
        data = HC.UniquePoints(points)
        @test HC.iscontained(data, [0.0im]) == false
    end

    @testset "Polynomials" begin
        @polyvar p q a b c x y z
        e = [p + 1]
        f = [a * b * c ]
        g = [x+y, y + z, x + z]
        @test expand(f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2]
        @test expand(e ∘ f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2 + 1]

        @test validate(f ∘ g) == true
        @test validate(f ∘ g, parameters=[z]) == true
        @test validate(f ∘ g, parameters=[c]) == false
        @test validate(g ∘ f) == false

        @test HC.nvariables(f ∘ g) == 3
        @test HC.nvariables(f ∘ g, parameters=[z]) == 2

        @polyvar x y z
        @test ishomogenous([x^2+y^2+x*y, x^5])
        @test ishomogenous([x^2+y^2+x*y, x^4+1]) == false

        #test weighted degree
        @test ishomogenous(x^3+x*y, [(x, 2), (y, 4)])
        @test homogenize(x+x*y, [(x, 2), (y, 4)], z) == x*z^4+x*y

        @test ishomogenous(homogenize([x^2+y^2+x*y, x^4+1]))
        @test ishomogenous([x^2+z^2 + y, x^4+z^4], [x,z]) == false
        @test ishomogenous(homogenize([x^2+z^2*y, x^4+z^4*y], [x,z]), [x,z]) == true

        @test ishomogenous(f ∘ g) == true
        h = [a * b * c  + p * c^3] ∘ [x+y, y + z, x + z]
        @test ishomogenous(h, parameters=[p]) == true
        h2 = [a * b * c  + p] ∘ [x+y, y + z, x + z]
        @test ishomogenous(h2, parameters=[p]) == false


        #homogenize
        @test homogenize(x^2+y+1, z) == x^2+y*z + z^2
        # This needs to be an array due to a compiler bug
        @test homogenize([x^2+y+p], z, parameters=[p]) == [x^2+y*z + z^2*p]

        @test homogenize([x^3+p, 1+y], z, parameters=[p]) == [x^3+p*z^3, z+y]

        h2 = [a * b * c  + p] ∘ [x+y, y + z, x + q]
        @test validate(homogenize(h2, parameters=[p, q]), parameters=[p, q])
        h3 = [a * b * c  + 1] ∘ [x+y, y + z, x + 1]
        @test validate(homogenize(h3))

        # Weylnorm
        @polyvar x y z
        f = 3.0x^2 + 2x*y - y^2
        g = (-2.5+2im) * x^2 - 3.0*x*y + 4y^2
        @test HC.weyldot(f, g) == 3.0 * conj(-2.5 + 2im) + 2.0 * (-3.0) / 2 + (-1.0) * 4.0
        @test HC.weyldot(f, f) == 9.0 + 4.0  / 2 + 1.0
        @test HC.weylnorm(f)^2 ≈ HC.weyldot(f,f)
        @test HC.weyldot([f, f], [g, g]) == 2 * HC.weyldot(f, g)
        @test HC.weylnorm([f, f]) == √HC.weyldot([f, f], [f, f])
    end

    @testset "Misc" begin
        A = rand(Complex{Float64}, 12, 12)
        b = rand(Complex{Float64}, 12)
        C, d = copy(A), copy(b)
        @test norm(HC.solve!(C, d) - A \ b) < 1e-10

        A = rand(Complex{Float64}, 15, 12)
        b = rand(Complex{Float64}, 15)
        C, d = copy(A), copy(b)
        HC.solve!(C, d)
        @test norm(d[1:12] - A \ b) < 1e-10
    end
end
