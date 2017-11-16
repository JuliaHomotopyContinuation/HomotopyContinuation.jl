export real_solutions, singular_solutions, isolated_solutions, solutions_at_infinity

"""
    real_solutions(R::HomotopyContinuation.Result)

Filters the real solution from the all the solutions in R.

"""
function real_solutions(ComplexResults::HomotopyContinuation.Result)

    RealResults = filter(S -> S.real_solution, ComplexResults)
    Result(RealResults)
end


"""
    singular_solutions(R::HomotopyContinuation.Result, tol)

Filters the singular solutions from the all the solutions in R.

"""

function singular_solutions(ComplexResults::HomotopyContinuation.Result)

    sol = filter(S -> S.returncode == :singular, ComplexResults)
    Result(sol)
end


"""
    isolated_solutions(R::HomotopyContinuation.Result, tol)

Filters the isolated solutions from the all the solutions in R.

"""

function isolated_solutions(ComplexResults::HomotopyContinuation.Result)

    sol = filter(S -> S.returncode == :isolated, ComplexResults)
    Result(sol)
end

"""
    solutions_at_infinity(R::HomotopyContinuation.Result, tol)

Filters the solutions at infinity from the all the solutions in R.

"""

function solutions_at_infinity(ComplexResults::HomotopyContinuation.Result)

    sol = filter(S -> S.returncode == :at_infinity, ComplexResults)
    Result(sol)
end

"""
    singular_at_infinity(R::HomotopyContinuation.Result, tol)

Filters the singular solutions at infinity from the all the solutions in R.

"""

function singular_at_infinity(ComplexResults::HomotopyContinuation.Result)

    sol = filter(S -> S.returncode == :singular_at_infinity, ComplexResults)
    Result(sol)
end
