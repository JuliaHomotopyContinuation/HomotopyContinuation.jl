@testset "Many parameters solver" begin
    # Setup
    @polyvar x y
    f = x^2 + y^2 - 1

    @polyvar a b c
    l = a * x + b * y + c
    F = [f, l]

    # Compute start solutions S₀ for given start parameters p₀
    p₀ = randn(ComplexF64, 3)
    S₀ = solutions(solve(subs(F, [a, b, c] => p₀)))
    # The parameters we are intersted in
    params = [rand(3) for i = 1:100]

    result1 = many_parameters_solve(
        F,
        S₀,
        p₀,
        params;
        parameters = [a, b, c],
        threading = true,
    )
    @test typeof(result1) == Vector{Tuple{Result{Vector{ComplexF64}},Vector{Float64}}}
    result1 = many_parameters_solve(
        F,
        S₀,
        p₀,
        params;
        parameters = [a, b, c],
        threading = false,
    )
    @test typeof(result1) == Vector{Tuple{Result{Vector{ComplexF64}},Vector{Float64}}}

    # Only keep real solutions
    result2 = many_parameters_solve(
        F,
        S₀,
        p₀,
        params;
        parameters = [a, b, c],
        transform_result = (r, p) -> real_solutions(r),
        threading = true,
    )
    @test typeof(result2) == Vector{Vector{Vector{Float64}}}
    @test !isempty(result2)

    # Now instead of an Array{Array{Array{Float64,1},1},1} we want to have an
    # Array{Array{Float64,1},1}
    result3 = many_parameters_solve(
        F,
        S₀,
        p₀,
        params;
        parameters = [a, b, c],
        transform_result = (r, p) -> real_solutions(r),
        flatten = true,
        threading = false,
    )
    @test typeof(result3) == Vector{Vector{Float64}}
    @test !isempty(result3)

    # The passed `params` do not directly need to be the target parameters.
    # Instead they can be some more concrete informations (e.g. an index)
    # and we can them by using the `transform_parameters` method
    result4 = many_parameters_solve(
        F,
        S₀,
        p₀,
        1:100;
        parameters = [a, b, c],
        transform_result = (r, p) -> (real_solutions(r), p),
        transform_parameters = _ -> rand(3),
    )
    @test typeof(result4) == Vector{Tuple{Vector{Vector{Float64}},Int64}}
end
