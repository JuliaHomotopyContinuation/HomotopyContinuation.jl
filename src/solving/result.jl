
"""
    Result(Vector{PathResult})

Constructs a summary of `PathResult`s.  A `Result` contains:
* `tracked`: Number of path that were tracked.
* `finite_solutions`: a vector containing all finite solutions.
* `solutions_at_infinity`: A vector containing all solutions at Infinity.
* `singular_solutions`: A vector containing all solution with windingnumber > 1.
* `failed_paths`: A vector containing all failed paths.
* `finite`: Length of `finite_solutions`
* `at_infinity`: Length of `solutions_at_infinity`
* `singular`:Length of `singular_solutions`
* `failed`: Length of `failed_paths`
"""

struct Result
    tracked::Int
    finite::Int
    at_infinity::Int
    singular::Int
    failed::Int

    finite_solutions::Vector{PathResult}
    solutions_at_infinity::Vector{PathResult}
    singular_solutions::Vector{PathResult}
    failed_paths::Vector{PathResult}
end

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

    Result(
    length(S),
    length(I),
    length(J),
    length(K),
    length(F),
    I, J, K, F
    )

end


function Base.show(io::IO, r::Result)

        println(io, "---------------------------------")
        println(io, "Paths tracked:  $(r.tracked) ($(r.failed) failed)")
        println(io, "# finite solutions:  $(r.finite)")
        println(io, "# solutions at Infinity:  $(r.at_infinity)")
        println(io, "# singular solutions (windingnumber > 1):  $(r.singular)")

        if r.finite > 0
            println(io, "---------------------------------")
            println(io, "Finite solutions:")
            for s in r.finite_solutions
                println(io, "$(round.(s.solution,4))")
            end
        end


        if r.at_infinity > 0
            println(io, "---------------------------------")
            println(io, "Solutions at infinity:")
            for s in r.solutions_at_infinity
                println(io, "$(round.(s.solution,4))")
            end
        end

        if r.singular > 0
            println(io, "---------------------------------")
            println(io, "Singular solutions (windingnumber > 1):")
            for s in r.singular_solutions
                println(io, "$(round.(s.solution,4))")
            end
        end
        println(io, "---------------------------------")
end


function Juno.render(i::Juno.Inline, r::Result)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths tracked → $(r.tracked) ($(r.failed) failed)"))]

        if r.finite > 0
                push!(t[:children],
                Juno.render(i, Text("Finite solutions:")),
                Juno.render(i, r.finite_solutions))
        end
        if r.at_infinity > 0
                push!(t[:children],
                Juno.render(i, Text("Solutions at ∞:")),
                Juno.render(i, r.solutions_at_infinity))
        end
        if r.singular > 0
                push!(t[:children],
                Juno.render(i, Text("Singular solutions:")),
                Juno.render(i, r.singular_solutions))
        end
        if r.failed > 0
                push!(t[:children],
                Juno.render(i, Text("Failed Paths:")),
                Juno.render(i, r.failed_paths))
        end

        return t
    end


"""
        real_solutions(Result)

    Filters the real solutions from `Result`. By default all solutions with absolute imaginary part of at most 1e-10 are returned. Optionally, one can pass a manually defined tolerance value `t` by `real_solutions(Result, t)`.
"""


function real_solutions(R::Result, tol::Float64)

        I = filter(
        s -> sqrt(maximum(abs2, imag.(s.solution))) < tol,
        R.finite_solutions)

        J = filter(
        s ->  sqrt(maximum(abs2, imag.(s.solution))) < tol,
        R.solutions_at_infinity)

        K = filter(
        s ->  sqrt(maximum(abs2, imag.(s.solution))) < tol,
        R.singular_solutions)

        F = filter(
        s ->  sqrt(maximum(abs2, imag.(s.solution))) < tol,
        R.failed_paths)

        Result(
        length(I) + length(J) + length(K) + length(F),
        length(I),
        length(J),
        length(K),
        length(F),
        I, J, K, F
        )
end
real_solutions(R::Solving.Result) = real_solutions(R, 1e-10)
