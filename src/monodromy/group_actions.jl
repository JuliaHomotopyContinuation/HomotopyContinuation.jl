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
