export ManyParametersSolver, many_parameters_solver, many_parameters_solve

"""
    ManyParametersSolver(solver)

A wrapper around a solver to prepare it for [`many_parameters_solve`](@ref).
"""
struct ManyParametersSolver{S<:Solver}
    solvers::Vector{S}
end
ManyParametersSolver(S::Solver) = ManyParametersSolver([S])


const many_parameters_solve_supported_keywords = [
    :transform_result,
    :transform_parameters,
    :flatten,
    :path_jumping_check,
    :threading,
]

"""
    many_parameters_solver( F,
         S::Vector{<:AbstractVector},
         start_params::AbstractVector{<:Number},
         params::AbstractVector; options...)

Construct a `ManyParametersSolver` by reusing the [`solver`](@ref) function.
"""
function many_parameters_solver(
    F,
    starts::Vector{<:AbstractVector},
    start_params::AbstractVector{<:Number},
    parameters::AbstractVector;
    kwargs...,
)
    ManyParametersSolver(solver(F, starts; generic_parameters = start_params, kwargs...))
end

"""
    many_parameters_solve(
        F,
        S::Vector{<:AbstractVector},
        start_params::AbstractVector{<:Number},
        params::AbstractVector;
        options...,
    )

Solve the system paramterized system `F` for many different parameters `parameters`.
If no further options are specified (resp. only those which can also be passed to [`solve`](@ref))
then the result of this funciton is similar to
```
map(params) do p
    res = solve(F, S; start_parameters = start_params, target_parameters = p, options...)
    (res, p)
end
```
Note if `threading = true` then the parallelization
is done on a parameter by parameter basis (and not as above on a path by path basis).
Note that order of the results is not guaranteed.
The `options` allow to apply transformations on the results

## Options
The function accepts most of the options which can be passed to [`solve`](@ref).
Furthermore, the following options can be passed.

* `flatten::Bool`: Flatten the output of `transform_result`. This is useful for example if
  `transform_result` returns a vector of solutions, and you only want a single vector of
  solutions as the result (instead of a vector of vector of solutions).
* `transform_result::Function`: A function taking two arguments, the `result` and the
  parameters `p`. By default this returns the tuple `(result, p)`.
* `transform_parameters::Function`: Transform a parameters values `p` before passing it to
  `target_parameters = ...`. By default this `identity`.
* `threading::Bool`: Enable multi-threading.

## Examples

We start with some setup. We intersect a circle with a many different lines.

```julia
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
```

Here are the some examples how to use the different options.
```julia-repl
julia> result1 = many_parameters_solve(F, S₀, p₀, params; parameters = [a, b, c]);

julia> typeof(result1)
Array{Tuple{Result{Array{Complex{Float64},1}},Array{Float64,1}},1}

julia> result1[1]
(Result{Array{Complex{Float64},1}} with 2 solutions
==================================================
• 2 non-singular solutions (2 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 217381
, [0.34822917261176745, 0.9745157077542976, 0.8920963081897844])

julia> # Only keep real solutions
       result2 = many_parameters_solve(
           F,
           S₀,
           p₀,
           params;
           parameters = [a, b, c],
           transform_result = (r, p) -> real_solutions(r),
       );

julia> typeof(result2)
Array{Array{Array{Float64,1},1},1}

julia> result2[1:3]
3-element Array{Array{Array{Float64,1},1},1}:
 [[-0.7673556041343146, -0.6412217844113383], [0.1872060592892163, -0.9823206662619913]]
 [[-0.5596149890742634, 0.8287527158347127], [0.5455922148592971, -0.8380507950505903]]
 [[-0.9836089512008734, 0.18031481114295062], [0.35645440225804775, -0.9343127201910812]]

julia> # Now instead of an Array{Array{Array{Float64,1},1},1} we want to have an
       # Array{Array{Float64,1},1}
       result3 = many_parameters_solve(
           F,
           S₀,
           p₀,
           params;
           parameters = [a, b, c],
           transform_result = (r, p) -> real_solutions(r),
           flatten = true
       );

julia> typeof(result3)
Array{Array{Float64,1},1}

julia> result3[1:3]
3-element Array{Array{Float64,1},1}:
 [-0.7673556041343146, -0.6412217844113383]
 [0.1872060592892163, -0.9823206662619913]
 [-0.5596149890742634, 0.8287527158347127]

julia> # The passed `params` do not directly need to be the target parameters.
       # Instead they can be some more concrete informations (e.g. an index)
       # and we can them by using the `transform_parameters` method
       result4 = many_parameters_solve(
           F,
           S₀,
           p₀,
           1:100;
           parameters = [a, b, c],
           transform_result = (r, p) -> (real_solutions(r), p),
           transform_parameters = _ -> rand(3)
       );

julia> typeof(result4)
Array{Tuple{Array{Array{Float64,1},1},Int64},1}

julia> result4[1:3]
3-element Array{Tuple{Array{Array{Float64,1},1},Int64},1}:
 ([[-0.9999999273168562, 0.0003812693043413411], [0.8944273254377234, -0.44721332662424196]], 1)
 ([[-0.9993469440132001, 0.03613427031887791], [0.8151981131530265, -0.5791822133938035]], 68)
 ([[-0.37073098818787914, 0.9287402943757951], [0.35361366239377473, -0.9353915638749696]], 69)
```
"""
function many_parameters_solve(
    F,
    starts::Vector{<:AbstractVector},
    start_params::AbstractVector{<:Number},
    parameters::AbstractVector;
    kwargs...,
)
    solve_kwargs, rest = splitkwargs(kwargs, many_parameters_solve_supported_keywords)
    solver = many_parameters_solver(F, starts, start_params, parameters; rest...)
    many_parameters_solve(solver, starts, start_params, parameters; solve_kwargs...)
end


function many_parameters_solve(
    MPS::ManyParametersSolver,
    starts::Vector{<:AbstractVector},
    start_params::AbstractVector{<:Number},
    parameters::AbstractVector;
    transform_result::Function = tuple,
    transform_parameters::Function = identity,
    flatten::Bool = false,
    path_jumping_check::Bool = true,
    threading::Bool = true,
)
    Threads.resize_nthreads!(MPS.solvers)

    start_parameters!.(MPS.solvers, Ref(start_params))

    first_result = solve(
        MPS.solvers[1],
        starts;
        target_parameters = transform_parameters(parameters[1]),
        threading = false,
        show_progress = false,
        path_jumping_check = path_jumping_check,
    )

    if flatten
        results = transform_result(first_result, parameters[1])
        if !(results isa AbstractArray)
            throw(ArgumentError("Cannot flatten arguments of type `$(typeof(results))`"))
        end
    else
        results = [transform_result(first_result, parameters[1])]
    end
    many_parameters_solve_barrier!(
        results,
        MPS,
        starts,
        parameters,
        transform_result,
        transform_parameters,
        Val(flatten),
    )
    results
end

function many_parameters_solve_barrier!(
    results,
    MPS::ManyParametersSolver,
    starts::Vector{<:AbstractVector},
    parameters::AbstractVector,
    transform_result::Function,
    transform_parameters::Function,
    ::Val{Flatten};
    path_jumping_check::Bool = false,
    threading::Bool = true,
) where {Flatten}
    N = length(parameters)

    if threading && Threads.nthreads() > 1
        thread_results = [similar(results, 0) for _ = 2:Threads.nthreads()]
        push!(thread_results, results)

        Threads.@threads for i = 2:N
            tid = Threads.threadid()
            rᵢ = transform_result(
                solve(
                    MPS.solvers[tid],
                    starts;
                    target_parameters = transform_parameters(parameters[i]),
                    path_jumping_check = path_jumping_check,
                    threading = false,
                    show_progress = false,
                ),
                parameters[i],
            )

            if Flatten
                append!(thread_results[tid], rᵢ)
            else
                push!(thread_results[tid], rᵢ)
            end
        end

        for i = 1:(Threads.nthreads()-1)
            append!(results, thread_results[i])
        end
    else
        for i = 2:N
            rᵢ = transform_result(
                solve(
                    MPS.solvers[1],
                    starts;
                    target_parameters = transform_parameters(parameters[i]),
                    path_jumping_check = path_jumping_check,
                    threading = false,
                    show_progress = false,
                ),
                parameters[i],
            )

            if Flatten
                append!(results, rᵢ)
            else
                push!(results, rᵢ)
            end
        end
    end

    results
end
