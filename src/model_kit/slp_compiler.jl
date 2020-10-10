function var_param_map_and_boundschecks(;
    nexpressions::Int,
    variables::Vector{Symbol},
    parameters::Vector{Symbol},
    t::Union{Symbol,Nothing} = nothing,
    # taylor = false,
    U::Bool = false,
    u::Bool = false,
)
    var_map = Dict{Symbol,Union{Symbol,Expr}}()
    for (i, v) in enumerate(variables)
        var_map[v] = :(x[$i])
    end
    for (i, v) in enumerate(parameters)
        var_map[v] = :(p[$i])
    end
    if t !== nothing
        var_map[t] = :(t)
    end
    m = length(variables)
    l = length(parameters)
    checks = Expr[]
    if U
        push!(checks, :(@boundscheck checkbounds(U, 1:$nexpressions, 1:$m)))
    end
    if u
        push!(checks, :(@boundscheck u === nothing || checkbounds(u, 1:$nexpressions)))
    end
    push!(checks, :(@boundscheck checkbounds(x, 1:$m)))
    push!(checks, :(@boundscheck p === nothing || checkbounds(p, 1:$l)))

    var_map, Expr(:block, checks...)
end

function unroll_pow(var, n)
    n == 0 && return :(one($var))
    n == 1 && return var
    n == 2 && return :(sqr($var))
    n == 3 && return :(@fastmath $var * sqr($var))
    n == 4 && return :(sqr(sqr($var)))
    n == 5 && return :(@fastmath $var * sqr(sqr($var)))

    # base to expansion shows up which power it is needed to compute
    d = digits(n, base = 2)
    x = :x
    exprs = [:(local $(Symbol(x, 1)) = $var)]
    for i = 2:length(d)
        push!(exprs, :(local $(Symbol(x, 1 << (i - 1))) = sqr($(Symbol(x, 1 << (i - 2))))))
    end
    prods = Symbol[]
    for (i, di) in enumerate(d)
        if !is_zero(di)
            push!(prods, Symbol(x, 1 << (i - 1)))
        end
    end
    if length(prods) > 1
        push!(exprs, :(@fastmath $(Expr(:call, :*, prods...))))
    end
    quote
        let
            $(exprs...)
        end
    end
end

function generate_evaluate_impl(
    list::InstructionList,
    out::Vector;
    variables::Vector{Symbol},
    t = nothing,
    parameters = Symbol[],
)
    var_param_map, boundschecks = var_param_map_and_boundschecks(
        nexpressions = length(out),
        variables = variables,
        parameters = parameters,
        t = t,
        u = true,
        U = false,
    )
    slp = to_expr(list; variables = var_param_map)
    assign = map(1:length(out), out) do i, a
        :(u[$i] = $(to_expr_arg(get(var_param_map, a, a))))
    end
    quote
        $boundschecks
        @inbounds begin
            $slp
            $(assign...)
        end
        u
    end
end

function generate_evaluate_and_jacobian_impl(
    list::InstructionList,
    out::Vector;
    variables::Vector{Symbol},
    t = nothing,
    parameters = Symbol[],
)
    grad_list, out, jac_out = ModelKit.gradient(list, out, variables)

    var_param_map, boundschecks = var_param_map_and_boundschecks(
        nexpressions = length(out),
        variables = variables,
        parameters = parameters,
        t = t,
        u = true,
        U = true,
    )
    slp = to_expr(grad_list; variables = var_param_map)
    jac_assign = Expr[]
    all_assigned = true
    for l = 1:length(variables), k = 1:length(out)
        a = jac_out[k, l]
        if a != 0
            push!(jac_assign, :(U[$k,$l] = $(to_expr_arg(get(var_param_map, a, a)))))
        else
            all_assigned = false
        end
    end

    quote
        $boundschecks
        $((!all_assigned ? (:(fill!(U, zero(eltype(U)))),) : ())...)
        @inbounds begin
            $slp
            if !(u isa Nothing)
                $((:(u[$i] = $(to_expr_arg(a))) for (i, a) in enumerate(out))...)
            end
            $(jac_assign...)
        end
        nothing
    end
end

function generate_taylor_impl(
    list::InstructionList,
    out::Vector;
    variables::Vector{Symbol},
    t = nothing,
    parameters = Symbol[],
)
    var_param_map, boundschecks = var_param_map_and_boundschecks(
        nexpressions = length(out),
        variables = variables,
        parameters = parameters,
        t = t,
        u = true,
        U = false,
    )
    slp = to_expr(list; variables = var_param_map)  do op, args
        Expr(:call, taylor_call_op(op), :(Val(K)), args...)
    end
    assign_taylor = map(1:length(out), out) do i, a
        :(u[$i] = $(to_expr_arg(get(var_param_map, a, a))))
    end
    assign_highest = map(1:length(out), out) do i, a
        :(u[$i] = last($(to_expr_arg(get(var_param_map, a, a)))))
    end
    if !isnothing(t)
        quote
            $boundschecks
            @inbounds let t = (t, one(t))
                $slp
                if u isa TaylorVector
                    $(assign_taylor...)
                else
                    $(assign_highest...)
                end
            end
            u
        end
    else
        quote
            $boundschecks
            @inbounds begin
                $slp
                if u isa TaylorVector
                    $(assign_taylor...)
                else
                    $(assign_highest...)
                end
            end
            u
        end
    end
end
