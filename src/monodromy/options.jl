struct Options{F1<:Function, F2<:Tuple, F3<:Function}
    target_solutions_count::Int
    tol::Float64
    timeout::Float64
    done_callback::F1
    group_actions::GroupActions{F2}
    parameter_sampler::F3
end

function Options(isrealsystem;
    target_solutions_count=error("target_solutions_count not provided"),
    tol::Float64=1e-6,
    timeout=float(typemax(Int)),
    done_callback=always_false,
    group_action=nothing,
    group_actions=GroupActions(group_action),
    parameter_sampler=independent_gaussian)

    if isrealsystem &&
       (group_actions isa GroupActions{Tuple{}}) # i.e. no group_actions provided
       actions = GroupActions(complex_conjugation)
    else
       actions = GroupActions(group_actions)
   end

    Options(target_solutions_count, tol, float(timeout), done_callback, actions, parameter_sampler)
end

always_false(x) = false
complex_conjugation(x) = (conj.(x),)


function independent_gaussian(p::SVector{N, T}) where {N, T}
    p + (@SVector randn(T, N))
end

function independent_zero_gaussian(p::SVector{N, T}) where {N, T}
    @SVector randn(T, N)
end
