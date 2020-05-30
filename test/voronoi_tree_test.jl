@testset "VoronoiTree" begin
    data = [randn(ComplexF64, 12) for i = 1:10_000]

    tree = HC.VoronoiTree(first(data), 1)
    for (i, d) in enumerate(data)
        insert!(tree, d, i)
    end

    @test length(tree) == length(data)
    @test sort!(collect(tree)) == 1:10_000
    @test all(i -> HC.search_in_radius(tree, data[i], 1e-12) == i, 1:length(data))
    @test all(enumerate(data)) do (i, d)
        HC.add!(tree, d, i, 1e-12) == i && HC.search_in_radius(tree, data[i], 1e-12) == i
    end
    @test all(enumerate(data)) do (i, d)
        HC.add!(tree, d, 0, 1e-12) == i
    end

    data = [randn(ComplexF64, 12) for i = 1:10_000]
    tree = HC.VoronoiTree(first(data), 1)
    for (i, d) in enumerate(data)
        insert!(tree, d, i)
    end
    d = data[9342] .+ 1e-5
    @test HC.search_in_radius(tree, d, 1e-4) == 9342
    @test HC.search_in_radius(tree, d, 1e-6) === nothing
end
