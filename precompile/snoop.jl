using HomotopyContinuation

@time begin
    @var x y
    affine_sqr =
        System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5])
    res = solve(affine_sqr; compile = false)
end
res = solve(affine_sqr; compile = false, start_system = :total_degree)
@var x y z
proj_square = System(
    [
        2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + (3 + 0im) * z^2,
        2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - (5 + 0im) * z^2,
    ]
)
solve(proj_square; compile = false)
solve(affine_sqr; compile = false, start_system = :total_degree)

@var x y
affine_ov = System(
    [
        (x^2 + y^2 + x * y - 3) * (x + 3),
        (x^2 + y^2 + x * y - 3) * (y - x + 2.1),
        2x + 5y - 3,
    ]
)
solve(affine_ov; compile = false)


function toric_ed(A)
    d, n = size(A)
    @var t[1:d] y[1:n] u[1:n]

    φ = map(j -> prod(i -> t[i]^A[i, j], 1:d), 1:n)
    Dφ = [differentiate(φ[i], t[j]) for i in 1:n, j in 1:d]

    return System([φ + y - u; Dφ' * y], parameters = u)
end

F = toric_ed([3 2 1 0; 0 1 2 3])
result = monodromy_solve(
    F,
    target_solutions_count = 21,
    max_loops_no_progress = 20,
    compile = false,
)

solve(
    F,
    solutions(result),
    start_parameters = parameters(result),
    target_parameters = randn(nparameters(F)),
    compile = false,
)
