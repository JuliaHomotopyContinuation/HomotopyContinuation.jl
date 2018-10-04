export monodromy_solve

import StaticArrays: SVector, @SVector
import LinearAlgebra


############
# Strategy
############

include("monodromy/strategy.jl")


################
# Group actions
################
"""
    GroupActions(actions::Function...)

Store a bunch of group actions `(f1, f2, f3, ...)`.
Each action has to return a tuple.
The actions are applied in the following sense
1) f1 is applied on the original solution `s`
2) f2 is applied on `s` and the results of 1
3) f3 is applied on `s` and the results of 1) and 2)
and so on

## Example
```julia-repl
julia> f1(s) = (s * s,);

julia> f2(s) = (2s, -s, 5s);

julia> f3(s) = (s + 1,);

julia> GroupActions(f1)(3)
(9,)

julia> GroupActions(f1,f2)(3)
(9, 18, -9, 45)

julia> GroupActions(f1,f2, f3)(3)
(9, 18, -9, 45, 10, 19, -8, 46)
```
"""
struct GroupActions{T<:Tuple}
    actions::T
end
GroupActions(::Nothing) = GroupActions(())
GroupActions(actions::GroupActions) = actions
GroupActions(actions::Function...) = GroupActions(actions)
function (group_action::GroupActions)(solution)
    apply_group_action(group_action.actions, solution)
end
apply_group_action(::Tuple{}, solution) = ()
# special case 1 function case
function apply_group_action(actions::Tuple{<:Function}, solution)
    (actions[1])(solution)
end
function apply_group_action(actions::Tuple, solution)
    f, rest = actions[1], Base.tail(actions)
    _apply_group_action(rest, f(solution))
end
function _apply_group_action(actions::Tuple, solutions::Tuple)
    foldl((acc, f) -> begin
        (acc..., foldl(flattentuple, map(f, acc); init=())...)
    end, actions; init=solutions)
end
flattentuple(a::Tuple, b::Tuple) = tuple(a..., b...)



##############
# Statistics
##############




struct MonodromyOptions{F1<:Function, F2<:Tuple}
    target_solutions_count::Int
    timeout::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
end

function MonodromyOptions(;
    target_solutions_count=error("target_solutions_count not provided"),
    timeout=float(typemax(Int)),
    done_callback=always_false,
    group_action=nothing,
    group_actions=GroupActions(group_action))

    MonodromyOptions(target_solutions_count, float(timeout), done_callback, GroupActions(group_actions))
end

always_false(x) = false
complex_conjugation(x) = (conj.(x),)




struct MonodromyStatistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    ntriangles::Int
end
MonodromyStatistics() = MonodromyStatistics(0, 0, 0)

function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike},
        p₀::AbstractVector{<:Number}, solution::Vector{<:Number}; kwargs...)

        monodromy_solve(F, p₀, [solution]; kwargs...)
end
function monodromy_solve(F::Vector{<:MP.AbstractPolynomialLike{T}},
        p₀::AbstractVector{<:Number},
        solutions::Vector{<:Vector{<:Number}};
        parameters=error("You need to provide `parameters=...` to monodromy"),
        strategy=Triangle(), options...) where {T}

    comp_solutions = map(v -> complex.(v), solutions)
    tracker = pathtracker(F, comp_solutions, parameters=parameters, p₁=p₀, p₀=p₀)
    # assemble strategy stuff
    nparams = length(p₀)
    p₀ = SVector{nparams}(p₀)
    strategy_params, strategy_cache = strategy_parameters_cache(strategy, tracker, p₀)
    monodromy_solve!(
        comp_solutions,
        tracker,
        p₀,
        strategy_params,
        strategy_cache,
        MonodromyOptions(;options...)
        )
end

function strategy_parameters_cache(strategy, tracker, p₀)
    parameters(strategy, p₀), cache(strategy, tracker)
end

function monodromy_solve!(
    solutions::Vector{<:Vector{<:Complex}},
    tracker::PathTracking.PathTracker,
    p₀::SVector,
    parameters::MonodromyStrategyParameters,
    strategy_cache::MonodromyStrategyCache,
    options::MonodromyOptions)

    t₀ = time_ns()
    k = 0

    # We prepopulate the solutions
    n = length(solutions)
    for i=1:n
        retcode = apply_group_actions_greedily!(nothing, solutions, solutions[i], options)
        if retcode == :done
            return :done
        end
    end


    while length(solutions) < options.target_solutions_count
        retcode = track_set!(solutions, tracker, p₀, parameters, strategy_cache, options)

        if retcode == :done
            break
        end

        dt = (time_ns() - t₀) * 1e-9
        if dt > options.timeout
            break
        end
        parameters = regenerate(parameters)
    end
    solutions
end

function track_set!(solutions, tracker, p₀, params::MonodromyStrategyParameters, strategy_cache, options::MonodromyOptions)
    S = copy(solutions)
    while !isempty(S)
        s₀ = pop!(S)
        s₁, retcode = loop(tracker, s₀, p₀, params, strategy_cache)
        if retcode == :success && !iscontained(solutions, s₁)
            push!(solutions, s₁)
            push!(S, s₁)
            if options.done_callback(s₁) || length(solutions) ≥ options.target_solutions_count
                return :done
            end

            retcode = apply_group_actions_greedily!(S, solutions, s₁, options)
            if retcode == :done
                return :done
            end
        end
    end
    :incomplete
end

function apply_group_actions_greedily!(S, solutions, s, options)
    for tᵢ in options.group_actions(s)
        if !iscontained(solutions, tᵢ)
            push!(solutions, tᵢ)
            if !(S isa Nothing)
                push!(S, tᵢ)
            end
            if options.done_callback(tᵢ) || length(solutions) ≥ options.target_solutions_count
                return :done
            end
            retcode = apply_group_actions_greedily!(S, solutions, tᵢ, options)
            if retcode == :done
                return :done
            end
        end
    end
    return :incomplete
end

function iscontained(solutions::Vector{T}, s_new; tol=1e-5) where {T<:AbstractVector}
    for s in solutions
        if LinearAlgebra.norm(s - s_new) < tol
        # if Distances.evaluate(Euclidean(), s, s_new) < tol
            return true
        end
    end
    false
end
