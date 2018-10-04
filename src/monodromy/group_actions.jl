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

(G::GroupActions)(solution) = apply_group_actions(G.actions, solution)

apply_group_actions(::Tuple{}, solution) = ()
# special case 1 function case
apply_group_actions(fs::Tuple{<:Function}, solution) = (fs[1])(solution)
# general case
function apply_group_actions(actions::Tuple, solution)
    f, rest = actions[1], Base.tail(actions)
    foldl(rest; init=f(solution)) do acc, f
        foldl(flatten, map(f, acc); init=acc)
    end
end
flatten(a::Tuple, b::Tuple) = (a..., b...)
