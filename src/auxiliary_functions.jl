export solutions, real_solutions, singular_solutions, isolated_solutions, solutions_at_infinity


"""
    solutions(r::Result; isolated=true, at_infnity=true, only_real=false, singular=true)

Filters the solutions which satisfy the constraints.
"""
function solutions(result::Result; isolated=true, at_infinity=true, only_real=false, singular=true)
    filter(result) do r
        if isolated && r.returncode == :isolated ||
           at_infinity && r.returncode == :at_infinity

            if only_real && singular
                return r.real_solution
            elseif only_real && !singular
                return r.real_solution && !r.singular
            elseif singular
                return true
            else
                return !r.singular
            end
        else
            return false
        end
    end
end


"""
    real_solutions(R::Result)

Filters the real solution from the all the solutions in R.

"""
function real_solutions(ComplexResults::Result)
    RealResults = filter(S -> S.real_solution, ComplexResults)
end


"""
    singular_solutions(R::Result)

Filters the singular solutions from the all the solutions in R.

"""

function singular_solutions(ComplexResults::Result)
    filter(S -> S.returncode == :singular, ComplexResults)
end


"""
    isolated_solutions(R::Result)

Filters the isolated solutions from the all the solutions in R.

"""

function isolated_solutions(ComplexResults::Result)
    filter(S -> S.returncode == :isolated, ComplexResults)
end

"""
    solutions_at_infinity(R::Result)

Filters the solutions at infinity from the all the solutions in R.

"""

function solutions_at_infinity(ComplexResults::Result)
    filter(S -> S.returncode == :at_infinity, ComplexResults)
end

"""
    singular_at_infinity(R::Result, tol)

Filters the singular solutions at infinity from the all the solutions in R.

"""

function singular_at_infinity(ComplexResults::Result)

    sol = filter(S -> S.returncode == :singular_at_infinity, ComplexResults)
    Result(sol)
end
