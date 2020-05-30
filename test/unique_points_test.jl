@testset "UniquePoints" begin
    data = [randn(ComplexF64, 12) for i = 1:10_000]
    up = HC.UniquePoints(first(data), 1)
    HC.add!.(up, data, 1:length(data), 1e-8)
    @test length(up) == 10_000
    HC.add!.(up, data, 1:length(data), 1e-8)


    # Test with group action
    x = randn(ComplexF64, 4)
    permutation1(x) = ([x[2]; x[1]; x[3]; x[4]],)
    permutation2(x) = ([x[1]; x[2]; x[4]; x[3]],)
    group_actions = GroupActions(permutation1, permutation1)
    X = [[x]; group_actions(x)]

    X = [v for v in GroupActions(permutation1, permutation2)(x)]

    # One group action
    data = HC.UniquePoints(first(X), 1, group_action = permutation1)
    HC.add!.(data, X, 1:length(X), 1e-5)
    @test length(data) == 2

    data =
        HC.UniquePoints(X[1], 1, group_actions = GroupActions(permutation1, permutation2))
    HC.add!.(data, X, 1:4, 1e-5)
    @test length(data) == 1

    # Group action and reality check
    x = randn(4)
    X = [im .* x, x, randn(ComplexF64, 4)]
    data = HC.UniquePoints(X[1], 1, group_action = x -> (im .* x, (-1) .* x, (-im) .* x))
    HC.add!.(data, X, 1:length(X), 1e-5)
    @test length(data) == 2
end
