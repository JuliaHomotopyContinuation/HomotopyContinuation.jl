
"""
    Result(Vector{PathResult})

Constructs a summary of `PathResult`s.  A `Result` contains:
* `PathResults` the vector of PathResults
* `tracked`: length of `PathResults`
* `finite`: Number of finite_solutions.
* `at_infinity`: Number of solutions at infinity.
* `singular`: Number of singular solutions.
* `failed`: Number of failed paths.
"""

struct Result
    PathResults::Vector{PathResult}
    tracked::Int
    finite::Int
    at_infinity::Int
    singular::Int
    failed::Int
end

function convert(Result, S::Vector{PathResult{T1, T2, T3}}) where {T1, T2, T3 <: Number}
    Result(
    S,
    length(S),
    count(s -> (s.returncode == :success) && s.windingnumber == 1, S),
    count(s -> (s.returncode == :at_infinity) && s.windingnumber == 1, S),
    count(s -> s.windingnumber > 1 && (s.returncode == :success || s.returncode == :at_infinity), S),
    count(s -> s.returncode == :path_failed, S)
    )
end
function convert(Result, S::Vector{PathResult})
    Result(
    S,
    length(S),
    count(s -> (s.returncode == :success) && s.windingnumber == 1, S),
    count(s -> (s.returncode == :at_infinity) && s.windingnumber == 1, S),
    count(s -> s.windingnumber > 1 && (s.returncode == :success || s.returncode == :at_infinity), S),
    count(s -> s.returncode == :path_failed, S)
    )
end


function finite_solutions(R::Result; compact = false)
    I = filter(
    s -> (s.returncode == :success) &&
    s.windingnumber == 1,
    R.PathResults)
    if !compact
        return convert(Result, I)
    else
        return map(s->s.solution, I)
    end
end

function solutions_at_infinity(R::Result; compact = false)
    J = filter(
    s -> (s.returncode == :at_infinity) &&
    s.windingnumber == 1,
    R.PathResults)
    if !compact
        return convert(Result, J)
    else
        return map(s->s.solution, J)
    end
end

function singular_solutions(R::Result; compact = false, tol = Inf)
    K = filter(
    s -> (s.windingnumber > 1 || s.condition_number > tol) &&
    (s.returncode == :success || s.returncode == :at_infinity),
    R.PathResults)
    if !compact
        return convert(Result, K)
    else
        return map(s->s.solution, K)
    end
end

function failed_paths(R::Result; compact = false)
    F = filter(
    s -> s.returncode == :path_failed,
    R.PathResults)
    if !compact
        return convert(Result, F)
    else
        return map(s->s.solution, F)
    end
end

function real_solutions(R::Result; compact = false, tol = 1e-10)
    real_sols = filter(
    s -> sqrt(maximum(abs2, imag.(s.solution))) < tol,
    R.PathResults)
    if !compact
        return convert(Result, real_sols)
    else
        return map(s->s.solution, real_sols)
    end
end



function Base.show(io::IO, r::Result)

        println(io, "-----------------------------------------------")
        println(io, "Paths:  $(r.tracked) ($(r.failed) failed)")
        println(io, "# finite solutions:  $(r.finite)")
        println(io, "# solutions at Infinity:  $(r.at_infinity)")
        println(io, "# singular solutions (windingnumber > 1):  $(r.singular)")

        if r.finite > 0
            println(io, "-----------------------------------------------")
            println(io, "Finite solutions:")
            for s in finite_solutions(r, compact = true)
                println(io, "$(round.(s,4))")
            end
        end


        if r.at_infinity > 0
            println(io, "-----------------------------------------------")
            println(io, "Solutions at infinity:")
            for s in solutions_at_infinity(r, compact = true)
                println(io, "$(round.(s,4))")
            end
        end

        if r.singular > 0
            println(io, "-----------------------------------------------")
            println(io, "Singular solutions (windingnumber > 1):")
            for s in singular_solutions(r, compact = true)
                println(io, "$(round.(s,4))")
            end
        end
        println(io, "-----------------------------------------------")
end


function Juno.render(i::Juno.Inline, r::Result)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths → $(r.tracked) ($(r.failed) failed)"))]

        if r.finite > 0
                push!(t[:children],
                Juno.render(i, Text("Finite solutions:")),
                Juno.render(i, round.(finite_solutions(r, compact = true), 4)))
        end
        if r.at_infinity > 0
                push!(t[:children],
                Juno.render(i, Text("Solutions at ∞:")),
                Juno.render(i, round.(solutions_at_infinity(r, compact = true), 4)))
        end
        if r.singular > 0
                push!(t[:children],
                Juno.render(i, Text("Singular solutions (windingnumber > 1):")),
                Juno.render(i, round.(singular_solutions(r, compact = true), 4)))
        end
        if r.failed > 0
                push!(t[:children],
                Juno.render(i, Text("Failed Paths:")),
                Juno.render(i, round.(failed_paths(r, compact = true), 4)))
        end
        return t
    end

#
#
