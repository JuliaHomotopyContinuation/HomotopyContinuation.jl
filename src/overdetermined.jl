export excess_solution_check, excess_solution_check!, OverdeterminedTracker
"""
    excess_solution_check!(path_result::PathResult,
                           F::RandomizedSystem,
                           newton_cache = NewtonCache(F.system))

Assigns to the [`PathResult`](@ref) `path_result` the `return_code` `:excess_solution` if
the `path_result` is a solution of the randomized system `F` but not of
the polynomial system underlying `F`.
This is performed by using Newton's method for non-singular solutions and comparing the
residuals of the solutions for singular solutions.
"""
function excess_solution_check!(
    path_result::PathResult,
    F::RandomizedSystem,
    newton_cache::NewtonCache = NewtonCache(F.system),
)
    is_success(path_result) || return path_result
    if is_nonsingular(path_result)
        acc = max(accuracy(path_result), eps())
        if path_result.extended_precision
            max_abs_norm_first_update = 100 * sqrt(acc)
        else
            max_abs_norm_first_update = 100acc
        end
        res = newton(
            F.system,
            solution(path_result),
            InfNorm(),
            newton_cache;
            extended_precision = path_result.extended_precision,
            abs_tol = 100 * acc,
            rel_tol = 100 * acc,
            max_abs_norm_first_update = max_abs_norm_first_update,
        )
        if !is_success(res)
            path_result.return_code = :excess_solution
            path_result.accuracy = res.accuracy
            path_result.residual = res.residual
        end
    else
        # for singular solutions compare residuals due to lack of something better right now
        evaluate!(newton_cache.r, F.system, solution(path_result))
        system_residual = LA.norm(newton_cache.r, InfNorm())
        if system_residual > 100 * residual(path_result)
            path_result.return_code = :excess_solution
            path_result.residual = system_residual
        end
    end

    path_result
end

struct ExcessSolutionCheck{S<:RandomizedSystem,M} <: Function
    system::S
    newton_cache::NewtonCache{M}
end
(check::ExcessSolutionCheck)(R::PathResult) =
    excess_solution_check!(R, check.system, check.newton_cache)

"""
    excess_solution_check(F::RandomizedSystem)

Returns a function `λ(::PathResult)` which performs the excess solution check.
The call `excess_solution_check(F)(path_result)` is identical to
`excess_solution_check!(F, path_result)`. See also [`excess_solution_check!`](@ref).
"""
excess_solution_check(F::RandomizedSystem) = ExcessSolutionCheck(F, NewtonCache(F.system))

"""
    OverdeterminedTracker(tracker::AbstractPathTracker, F::RandomizedSystem)

Wraps the given [`AbstractPathTracker`](@ref) `tracker` to apply
[`excess_solution_check`](@ref) for the given randomized system `F` on each path result.
"""
struct OverdeterminedTracker{T<:AbstractPathTracker,E<:ExcessSolutionCheck} <:
       AbstractPathTracker
    tracker::T
    excess_solution_check::E
end
OverdeterminedTracker(T::AbstractPathTracker, F::RandomizedSystem) =
    OverdeterminedTracker(T, excess_solution_check(F))

track(T::OverdeterminedTracker, x; kwargs...) =
    T.excess_solution_check(track(T.tracker, x; kwargs...))
