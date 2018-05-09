
"""
    Result(Vector{PathResult})

Constructs a summary of `PathResult`s.  A `Result` contains:
* `solutions` the vector of PathResults
* `tracked`: length of `solutions`
* `finite_solutions(; compact = false)`: returns all finite solutions. `compact = true` returns only the computed roots with no additional information.
* `solutions_at_infinity(; compact = false)`: returns all solutions at Infinity.
* `singular_solutions(; compact = false)`: returns all solution with windingnumber > 1.
* `real_solutions(; compact = false, tol = 1e-10)`: returns all real solutions, that is all points with imaginary part at most `tol`.
* `failed_paths()`: returns all failed paths.
* `finite`: Length of `finite_solutions()`
* `at_infinity`: Length of `solutions_at_infinity()`
* `singular`:Length of `singular_solutions()`
* `failed`: Length of `failed_paths()`
"""

struct Result
    tracked::Int
    finite::Int
    at_infinity::Int
    singular::Int
    failed::Int

    solutions::Vector{PathResult}

    finite_solutions::Function
    solutions_at_infinity::Function
    singular_solutions::Function
    real_solutions::Function
    failed_paths::Function

    function Result(S::Vector{PathResult{T1, T2, T3}}) where {T1, T2, T3 <: Number}
        I = filter(
        s -> (s.returncode == :success) &&
        s.windingnumber == 1,
        S)

        J = filter(
        s -> (s.returncode == :at_infinity) &&
        s.windingnumber == 1,
        S)

        K = filter(
        s -> s.windingnumber > 1 &&
        (s.returncode == :success || s.returncode == :at_infinity),
        S)

        F = filter(
        s -> s.returncode == :path_failed,
        S)

        function finite_solutions(; compact = false)
            if !compact
                return Result(I)
            else
                return map(s->s.solution, I)
            end
        end

        function solutions_at_infinity(; compact = false)
            if !compact
                return Result(J)
            else
                return map(s->s.solution, J)
            end
        end

        function singular_solutions(; compact = false)
            if !compact
                return Result(K)
            else
                return map(s->s.solution, K)
            end
        end

        function failed_paths(; compact = false)
            if !compact
                return Result(F)
            else
                return map(s->s.solution, F)
            end
        end

        function real_solutions(; compact = false, tol = 1e-10)
            real_sols = filter(
            s -> sqrt(maximum(abs2, imag.(s.solution))) < tol,
            R)
            if !compact
                return real_sols
            else
                return map(s->s.solution, real_sols)
            end
        end

        new(
        length(S),
        length(I),
        length(J),
        length(K),
        length(F),
        S,
        finite_solutions,
        solutions_at_infinity,
        singular_solutions,
        real_solutions,
        failed_paths
        )
    end
end




function Base.show(io::IO, r::Result)

        println(io, "-----------------------------------------------")
        println(io, "Paths tracked:  $(r.tracked) ($(r.failed) failed)")
        println(io, "# finite solutions:  $(r.finite)")
        println(io, "# solutions at Infinity:  $(r.at_infinity)")
        println(io, "# singular solutions (windingnumber > 1):  $(r.singular)")

        if r.finite > 0
            println(io, "-----------------------------------------------")
            println(io, "Finite solutions:")
            for s in r.finite_solutions(; compact = true)
                println(io, "$(round.(s,4))")
            end
        end


        if r.at_infinity > 0
            println(io, "-----------------------------------------------")
            println(io, "Solutions at infinity:")
            for s in r.solutions_at_infinity(; compact = true)
                println(io, "$(round.(s,4))")
            end
        end

        if r.singular > 0
            println(io, "-----------------------------------------------")
            println(io, "Singular solutions (windingnumber > 1):")
            for s in r.singular_solutions(; compact = true)
                println(io, "$(round.(s,4))")
            end
        end
        println(io, "-----------------------------------------------")
end


function Juno.render(i::Juno.Inline, r::Result)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths tracked → $(r.tracked) ($(r.failed) failed)"))]

        if r.finite > 0
                push!(t[:children],
                Juno.render(i, Text("Finite solutions:")),
                Juno.render(i, r.finite_solutions(; compact = true)))
        end
        if r.at_infinity > 0
                push!(t[:children],
                Juno.render(i, Text("Solutions at ∞:")),
                Juno.render(i, r.solutions_at_infinity(; compact = true)))
        end
        if r.singular > 0
                push!(t[:children],
                Juno.render(i, Text("Singular solutions:")),
                Juno.render(i, r.singular_solutions(; compact = true)))
        end
        if r.failed > 0
                push!(t[:children],
                Juno.render(i, Text("Failed Paths:")),
                Juno.render(i, r.failed_paths(; compact = true)))
        end
        return t
    end

#
#
