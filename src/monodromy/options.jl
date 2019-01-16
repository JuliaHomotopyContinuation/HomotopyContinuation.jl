const options_allowed_keywords = [:tol, :done_callback,
    :group_action,:group_actions, :group_action_on_all_nodes,
    :parameter_sampler, :target_solutions_count, :timeout,
    :minimal_number_of_solutions, :maximal_number_of_iterations_without_progress]

struct Options{F1<:Function, F2<:Tuple, F3<:Function}
    tol::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
    group_action_on_all_nodes::Bool
    parameter_sampler::F3
    # stopping heuristic
    target_solutions_count::Int
    timeout::Float64
    minimal_number_of_solutions::Int
    maximal_number_of_iterations_without_progress::Int
end

function Options(isrealsystem;
    tol::Float64=1e-5,
    done_callback=always_false,
    group_action=nothing,
    group_actions=GroupActions(group_action),
    group_action_on_all_nodes=false,
    parameter_sampler=independent_normal,
    # stopping heuristic
    target_solutions_count=nothing,
    timeout=float(typemax(Int)),
    minimal_number_of_solutions::Int=default_minimal_number_of_solutions(target_solutions_count),
    maximal_number_of_iterations_without_progress::Int=10)

    if isrealsystem &&
       (group_actions isa GroupActions{Tuple{}}) # i.e. no group_actions provided
       actions = GroupActions(complex_conjugation)
    else
       actions = GroupActions(group_actions)
   end

    Options(tol, done_callback, actions,
        group_action_on_all_nodes, parameter_sampler,
        target_solutions_count == nothing ? typemax(Int) : target_solutions_count,
        float(timeout),
        minimal_number_of_solutions,
        maximal_number_of_iterations_without_progress)
end

default_minimal_number_of_solutions(::Nothing) = 2
function default_minimal_number_of_solutions(target_solutions_count::Int)
    div(target_solutions_count, 2)
end

always_false(x, sols) = false

"""
    complex_conjugation(x)

A group action which returns the elementwise complex conjugated solutions.
"""
complex_conjugation(x) = (conj.(x),)
has_group_actions(options::Options) = !(options.group_actions isa GroupActions{Tuple{}})

function independent_gaussian(p::SVector{N, T}) where {N, T}
    p + (@SVector randn(T, N))
end

function independent_normal(p::SVector{N, T}) where {N, T}
    @SVector randn(T, N)
end
