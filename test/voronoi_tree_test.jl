@testset "VoronoiTree" begin
    data = [randn(ComplexF64, 12) for i = 1:10_000]

    tree = HC.VoronoiTree(first(data), 1)
    for (i, d) in enumerate(data)
        insert!(tree, d, i)
    end

    @test length(tree) == length(data)
    @test sort!(collect(tree)) == 1:10_000
    @test all(i -> HC.search_in_radius(tree, data[i], 1e-12) == i, 1:length(data))
    let tree = HC.VoronoiTree(first(data), 1)
        @test all(enumerate(data)) do (i, d)
            HC.add!(tree, d, i, 1e-12) == (i, true) &&
                HC.search_in_radius(tree, d, 1e-12) == i
        end
    end

    data = [randn(ComplexF64, 12) for i = 1:10_000]
    tree = HC.VoronoiTree(first(data), 1)
    for (i, d) in enumerate(data)
        insert!(tree, d, i)
    end
    d = data[9342] .+ 1e-5
    @test HC.search_in_radius(tree, d, 1e-4) == 9342
    @test HC.search_in_radius(tree, d, 1e-6) === nothing

    # Test many points with nearly indentical distance to the inserted point
    p = shuffle!([[cis(k / 100 * 2π)] for k = 0:99])
    tree = HC.VoronoiTree(p[1], 1)
    for (i, pi) in enumerate(p)
        insert!(tree, pi, i)
    end
    @test isnothing(HC.search_in_radius(tree, [0.0], 1e-5))
    insert!(tree, [1e-5], 101)
    @test isnothing(HC.search_in_radius(tree, [0.0], 1e-6))
    @test HC.search_in_radius(tree, [0.0], 1e-4) == 101
end
