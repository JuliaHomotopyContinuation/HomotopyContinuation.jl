using HomotopyContinuation, Test

d = 2
@var x y a[1:6]
F = System(
    [
        (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
        (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
    ];
    parameters = a,
)
res = solve(
    F;
    iterator_only = true,
    target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32],
)


@test startswith(sprint(show, res), "ResultIterator")
@test seed(res) isa UInt32

@test length(path_results(res)) == ntracked(res) == 7
@test length(collect(results(res))) == nresults(res) == 3
#code length(ri::Base.Generators...blah blah but we know it has to do with a ResultIterator) to throw an error and advice. 
@test typeof(solutions(res)) <: Base.Generator

real_iter = real_solutions(res)
@test typeof(real_iter) <: Base.Generator
@test collect(real_iter) isa Vector
@test count(x -> true, real_iter) == nreal(res) == 1
real_pos = imap(x -> all(x .>  0), real_iter) |> collect
@test real_pos[1]

@test nnonsingular(res) == 3
@test isempty(singular(res))
@test nsingular(res) == 0
@test nexcess_solutions(res) == 0