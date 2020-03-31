## Variable construction

Variable(name::Union{Symbol,AbstractString}, indices::Int...) =
    Variable("$(name)$(join(map_subscripts.(indices), "₋"))")

const SUBSCRIPTS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
const SUBSCRIPT_MAP = Dict([first(string(i)) => SUBSCRIPTS[i+1] for i = 0:9])
const SUBSCRIPT_TO_INT_MAP = Dict([SUBSCRIPTS[i+1] => i for i = 0:9])
map_subscripts(index) = join(SUBSCRIPT_MAP[c] for c in string(index))

function sort_key(v::Variable)
    name = string(ModelKit.name(v))
    sub_start = findnext(c -> c in ModelKit.SUBSCRIPTS, name, 1)
    if isnothing(sub_start)
        return name, Int[]
    end
    var_base = name[1:sub_start-1]
    str_indices = split(name[sub_start:end], "₋")
    indices = map(str_indices) do str_index
        digits = map(s -> ModelKit.SUBSCRIPT_TO_INT_MAP[s], collect(str_index))
        n = length(digits)
        sum([digits[n-k+1] * 10^(k - 1) for k = 1:n])
    end
    var_base, reverse(indices)
end
Base.isless(a::Variable, b::Variable) = isless(sort_key(a), sort_key(b))
Base.sort!(vs::AbstractVector{Variable}; kwargs...) =
    permute!(vs, sortperm(sort_key.(vs); kwargs...))

Symbol(v::Variable) = name(v)


"""
    @var variable1 variable2 ...

Declare variables with the given names and automatically create the variable bindings.
The macro supports indexing notation to create `Array`s of variables.

## Examples

```julia
julia> @var a b x[1:2] y[1:2,1:3]
(a, b, Variable[x₁, x₂], Variable[y₁₋₁ y₁₋₂ y₁₋₃; y₂₋₁ y₂₋₂ y₂₋₃])

julia> a
a

julia> b
b

julia> x
2-element Array{Variable,1}:
 x₁
 x₂

julia> y
2×3 Array{Variable,2}:
 y₁₋₁  y₁₋₂  y₁₋₃
 y₂₋₁  y₂₋₂  y₂₋₃
```
"""
macro var(args...)
    vars, exprs = buildvars(args; unique = false)
    :($(foldl((x, y) -> :($x; $y), exprs, init = :()));
    $(Expr(:tuple, esc.(vars)...)))
end

"""
    @unique_var(args...)

Declare variables and automatically create the variable bindings to the given names.
This will change the names of the variables to ensure uniqueness.

## Examples

```julia
julia> @unique_var a b
(##a#591, ##b#592)

julia> a
##a#591

julia> b
##b#592
```
"""
macro unique_var(args...)
    vars, exprs = buildvars(args; unique = true)
    :($(foldl((x, y) -> :($x; $y), exprs, init = :()));
    $(Expr(:tuple, esc.(vars)...)))
end

function var_array(prefix, indices...)
    map(i -> Variable(prefix, i...), Iterators.product(indices...))
end

function buildvar(var; unique::Bool = false)
    if isa(var, Symbol)
        varname = unique ? gensym(var) : var
        var, :($(esc(var)) = Variable($"$varname"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) ||
        error("Expected $var to be of the form varname[idxset]")
        (2 ≤ length(var.args)) || error("Expected $var to have at least one index set")
        varname = var.args[1]
        prefix = unique ? string(gensym(varname)) : string(varname)
        varname, :($(esc(varname)) = var_array($prefix, $(esc.(var.args[2:end])...)))
    end
end

function buildvars(args; unique::Bool = false)
    vars = Symbol[]
    exprs = []
    for arg in args
        if arg isa Expr && arg.head == :tuple
            for _arg in arg.args
                var, expr = buildvar(_arg; unique = unique)
                push!(vars, var)
                push!(exprs, expr)
            end
        else
            var, expr = buildvar(arg; unique = unique)
            push!(vars, var)
            push!(exprs, expr)
        end
    end
    vars, exprs
end

"""
    unique_variable(var::Symbol, variables, parameters)

Create a variable with a unique name which doesn't clash with `variables` or
`parameters`.
If `var` is not possible the names `var##k` for `k=0,1,...` are tried until
one is possible,
"""
function unique_variable(var::Symbol, vars::Vector{Variable}, params::Vector{Variable})
    v = Variable(var)
    k = 0
    while (v in vars || v in params)
        v = Variable("$var##$k")
        k += 1
    end
    v
end


Base.adjoint(expr::Basic) = expr
Base.conj(expr::Basic) = expr
Base.transpose(expr::Basic) = expr
Base.broadcastable(v::Basic) = Ref(v)

"""
    variables(expr::Expression, parameters = Variable[])
    variables(exprs::AbstractVector{<:Expression}, parameters = Variable[])

Obtain all variables used in the given expression as an `Set`.
"""
function variables(exprs::AbstractVector{<:Basic})
    S = Set{Variable}()
    for expr in exprs
        union!(S, variables(expr))
    end
    sort!(collect(S))
end
function variables(exprs::Union{Basic,AbstractVector{<:Basic}}, params)
    setdiff!(variables(exprs), params)
end

"""
    nvariables(expr::Expression, parameters = Variable[])
    nvariables(exprs::AbstractVector{<:Expression}, parameters = Variable[])

Obtain the number of variables used in the given expression.
"""
nvariables(E::Union{Expression,AbstractVector{<:Expression}}, p = Variable[]) =
    length(variables(E, p))


"""
    subs(expr::Expression, subs::Pair{Variable,<:Expression}...)
    subs(expr::Expression, subs::Pair{AbstractArray{<:Variable},AbstractArray{<:Expression}}...)

Substitute into the given expression.

## Example

```
@var x y z

julia> subs(x^2, x => y)
y ^ 2

julia> subs(x * y, [x,y] => [x+2,y+2])
(x + 2) * (y + 2)
```
"""
subs(ex::Basic, args...) = subs(ex, ExpressionMap(), args...)
subs(ex::Basic, D::Dict) = subs(ex, ExpressionMap(D))
function subs(
    ex::Basic,
    D::ExpressionMap,
    (xs, ys)::Pair{<:AbstractArray{<:Basic},<:AbstractArray},
    args...,
)
    size(xs) == size(ys) ||
    throw(ArgumentError("Substitution arguments don't have the same size."))
    for (x, y) in zip(xs, ys)
        D[x] = y
    end
    subs(ex, D, args...)
end
function subs(ex::Basic, D::ExpressionMap, (x, y), args...)
    D[Expression(x)] = Expression(y)
    subs(ex, D, args...)
end
subs(exs::AbstractArray{<:Basic}, args...) = map(ex -> subs(ex, args...), exs)

"""
    evaluate(expr::Expression, subs...)
    evaluate(expr::AbstractArray{<:Expression}, subs...)

Evaluate the given expression.

## Example

```
@var x y

julia> evaluate(x^2, x => 2)
4

julia> evaluate(x * y, [x,y] => [2, 3])
6
"""
function evaluate(expr::AbstractArray{<:Basic}, args...)
    out = map(to_number, subs(expr, args...))
    if eltype(out) in (Number, Any)
        return to_smallest_eltype(out)
    else
        out
    end
end
evaluate(expr::Basic, args...) = to_number(subs(expr, args...))
(f::Union{Basic,AbstractArray{<:Basic}})(args...) = evaluate(f, args...)

"""
    evaluate!(u, expr::AbstractArray{<:Expression}, subs...)

Inplace for of [`evaluate`](@ref).
"""
function evaluate!(u::AbstractArray, expr::AbstractArray{<:Basic}, args...)
    map!(to_number, u, subs(expr, args...))
end

"""
    differentiate(expr::Expression, var::Variable)
    differentiate(expr::Expression, var::Vector{Variable})
    differentiate(expr::::Vector{<:Expression}, var::Variable, k = 1)
    differentiate(expr::Vector{<:Expression}, var::Vector{Variable})

Compute the derivative of `expr` with respect to the given variable `var`.
"""
function differentiate(expr::Basic, vars::AbstractVector{Variable})
    [differentiate(expr, v) for v in vars]
end
function differentiate(exprs::AbstractVector{<:Basic}, var::Variable, k = 1)
    [differentiate(e, var, k) for e in exprs]
end
function differentiate(exprs::AbstractVector{<:Basic}, vars::AbstractVector{Variable})
    [differentiate(e, v) for e in exprs, v in vars]
end

"""
    monomials(vars::Vector{<:Variable}, d; homogeneous::Bool = false)

Create all monomials of a given degree.

```
julia> @var x y
(x, y)

julia> monomials([x,y], 2)
6-element Array{Expression,1}:
   1
   x
   y
 x^2
 x*y
 y^2

julia> monomials([x,y], 2; homogeneous = true)
3-element Array{Operation,1}:
 x ^ 2
 x * y
 y ^ 2
 ```
"""
function monomials(vars::AbstractVector{Variable}, d::Integer; homogeneous::Bool = false)
    n = length(vars)
    if homogeneous
        pred = x -> sum(x) == d
    else
        pred = x -> sum(x) ≤ d
    end
    exps = collect(Iterators.filter(pred, Iterators.product(Iterators.repeated(0:d, n)...)))
    sort!(exps, lt = td_order)
    map(exps) do exp
        prod(i -> vars[i]^exp[i], 1:n)
    end
end

function rand_poly(vars::AbstractVector{Variable}, d::Integer; kwargs...)
    rand_poly(ComplexF64, vars, d; kwargs...)
end
function rand_poly(T, vars::AbstractVector{Variable}, d::Integer; homogeneous::Bool = false)
    M = monomials(vars, d; homogeneous = homogeneous)
    sum(randn(T, length(M)) .* M)
end

function td_order(x, y)
    sx = sum(x)
    sy = sum(y)
    sx == sy ? x > y : sx < sy
end

"""
    expand(e::Expression)

Expand a given expression.

```julia
julia> @var x y
(x, y)

julia> expand((x + y) ^ 2)
2*x*y + x^2 + y^2
```
"""
expand(e::Basic) = symengine_expand(e)

function to_dict(expr::Expression, vars::AbstractVector{Variable})
    mul_args, pow_args = ExprVec(), ExprVec()
    dict = Dict{Vector{Int},Expression}()

    if class(expr) == :Add
        for op in args(expr)
            to_dict_op!(dict, op, vars, mul_args, pow_args)
        end
    else
        to_dict_op!(dict, expr, vars, mul_args, pow_args)
    end

    dict
end

function to_dict_op!(dict, op, vars, mul_args, pow_args)
    cls = class(op)
    d = zeros(Int, length(vars))
    coeff = Expression(1)
    if cls == :Mul
        op_args = args!(mul_args, op)
        for arg in op_args
            arg_cls = class(arg)
            is_coeff = true
            if arg_cls == :Symbol
                for (i, v) in enumerate(vars)
                    if v == arg
                        d[i] = 1
                        is_coeff = false
                        break
                    end
                end
            elseif arg_cls == :Pow
                vec = args!(pow_args, arg)
                x = vec[1]
                for (i, v) in enumerate(vars)
                    if x == v
                        d[i] = convert(Int, vec[2])
                        is_coeff = false
                        break
                    end
                end
            end
            if is_coeff
                mul!(coeff, coeff, arg)
            end
        end
    elseif cls == :Symbol
        is_in_vars = false
        for (i, v) in enumerate(vars)
            if v == op
                d[i] = 1
                is_in_vars = true
                break
            end
        end
        if !is_in_vars
            coeff = copy(op)
        end
    elseif cls == :Pow
        vec = args!(pow_args, op)
        x = vec[1]
        for (i, v) in enumerate(vars)
            if x == v
                d[i] = convert(Int, vec[2])
                break
            end
        end
    elseif cls == :Div
        div!(coeff, coeff, op)
    elseif cls ∈ NUMBER_TYPES
        coeff = copy(op)
    end

    if haskey(dict, d)
        add!(dict[d], dict[d], coeff)
    else
        dict[d] = coeff
    end
    dict
end

"""
    exponents_coefficients(f::Expression, vars::AbstractVector{Variable}; expanded = false)

Return a matrix `M` containing the exponents for all occuring terms
(one term per column) and a vector `c` containing the corresponding coefficients.
Expands the given expression `f` unless `expanded = true`.
"""
function exponents_coefficients(
    f::Expression,
    vars::AbstractVector{Variable};
    expanded::Bool = false,
)
    expanded || (f = expand(f))
    D = to_dict(f, vars)
    # make exponents to matrix
    m = length(D)
    n = length(vars)
    M = zeros(Int, n, m)
    for (j, E) in enumerate(keys(D))
        for (i, dᵢ) in enumerate(E)
            M[i, j] = dᵢ
        end
    end
    coeffs = collect(values(D))
    M, coeffs
end

"""
    coefficients(f::Expression, vars::AbstractVector{Variable}; expanded = false)

Return all coefficients of the given polynomial `f` for the given variables `vars`.
"""
function coefficients(f::Expression, vars::AbstractVector{Variable}; expanded::Bool = false)
    expanded || (f = expand(f))
    collect(values(to_dict(f, vars)))
end

"""
    degree(f::AbstractVector{Expression}, vars = variables(f); expanded = false)

Compute the degrees of the expressions `f` in `vars`.
Unless `expanded` is `true` the expressions are first expanded.
"""
function degree(
    f::AbstractVector{Expression},
    vars::AbstractVector{Variable} = variables(f);
    expanded::Bool = false,
)
    if !expanded
        f = ModelKit.expand.(f)
    end
    dicts = ModelKit.to_dict.(f, Ref(vars))
    maximum.(sum, keys.(dicts))
end

"""
    degree(f::Expression, vars = variables(f); expanded = false)

Compute the degree of the expression `f`  in `vars`.
Unless `expanded` is `true` the expression is first expanded.
"""
function degree(
    f::Expression,
    vars::AbstractVector{Variable} = variables(f);
    expanded::Bool = false,
)
    if !expanded
        f = ModelKit.expand(f)
    end
    dicts = ModelKit.to_dict(f, vars)
    maximum(sum, keys(dicts))
end

function LinearAlgebra.det(M::AbstractMatrix{<:Union{Variable,Expression}})
    m = size(M)[1]
    if m > 2
        return sum(
            (-1)^(i - 1) * M[i, 1] * LinearAlgebra.det(M[1:end.!=i, 2:end]) for i = 1:m
        )
    else
        return M[1, 1] * M[2, 2] - M[2, 1] * M[1, 2]
    end
end

################
# Optimization #
################

"""
    horner(f::Expression, vars = variables(f))

Rewrite `f` using a multi-variate horner schema.

### Example
```julia
julia> @var u v c[1:3]
(u, v, Variable[c₁, c₂, c₃])

julia> f = c[1] + c[2] * v + c[3] * u^2 * v^2 + c[3]u^3 * v
c₁ + v*c₂ + u^2*v^2*c₃ + u^3*v*c₃

julia> ModelKit.horner(f)
c₁ + v*(c₂ + u^3*c₃ + u^2*v*c₃)
```
"""
function horner(f::Expression, vars = variables(f))
    M, coeffs = ModelKit.exponents_coefficients(f, vars)
    multivariate_horner(M, coeffs, vars)
end

function horner(coeffs::AbstractVector{Expression}, var::Variable)
    d = length(coeffs) - 1
    h = copy(coeffs[d+1])
    @inbounds for k = d:-1:1
        ModelKit.mul!(h, h, var)
        ModelKit.add!(h, h, coeffs[k])
    end
    h
end

function multivariate_horner(
    M::AbstractMatrix{Int},
    coeffs::AbstractVector{Expression},
    vars::AbstractVector{Variable},
)
    n, m = size(M)
    if m == 1
        fast_c = coeffs[1]
        for (i, var) in enumerate(vars)
            if M[i, 1] != 0
                if M[i, 1] == 1
                    fast_c *= var
                else
                    fast_c *= var^M[i, 1]
                end
            end
        end
        return fast_c::Expression
    elseif n == 1
        d = maximum(M)
        uni_coeff_indices = Union{Nothing,Int}[nothing for _ = 0:d]
        for j = 1:size(M, 2)
            uni_coeff_indices[M[1, j]+1] = j
        end

        var_coeffs = map(uni_coeff_indices) do ind
            if isnothing(ind)
                return zero(vars[1])
            else
                coeffs[ind]
            end
        end

        horner(var_coeffs, vars[1])
    else
        counts = zeros(Int, n)
        @inbounds for j = 1:size(M, 2), i = 1:size(M, 1)
            counts[i] += M[i, j] > 0
        end

        _, var_ind = findmax(counts)
        var = vars[var_ind]

        # compute degree in vars[var_ind]
        d = 0
        for j = 1:size(M, 2)
            d = max(d, M[var_ind, j])
        end

        # find coefficients
        coeff_indices = [Int[] for _ = 0:d]
        for j = 1:size(M, 2)
            push!(coeff_indices[M[var_ind, j]+1], j)
        end

        reduced_vars_ind = Int[]
        for i = 1:n
            if counts[i] > 0 && i != var_ind
                push!(reduced_vars_ind, i)
            end
        end
        reduced_vars = view(vars, reduced_vars_ind)

        var_coeffs = map(coeff_indices) do ind
            isempty(ind) && return zero(var)
            if length(ind) == 1
                j = first(ind)
                c = coeffs[j]
                for i in reduced_vars_ind
                    if M[i, j] != 0
                        if M[i, j] == 1
                            c *= vars[i]
                        else
                            c *= vars[i]^M[i, j]
                        end
                    end
                end
                return c::Expression
            end
            M_ind = view(M, reduced_vars_ind, ind)
            coeffs_ind = view(coeffs, ind)
            multivariate_horner(M_ind, coeffs_ind, reduced_vars)::Expression
        end
        horner(var_coeffs, var)
    end
end



#########################
## System and Homotopy ##
#########################


############
## System ##
############

function check_vars_params(f, vars, params)
    vars_params = params === nothing ? vars : [vars; params]
    Δ = setdiff(variables(f), vars_params)
    isempty(Δ) || throw(ArgumentError(
        "Not all variables or parameters of the system are given. Missing: " *
        join(Δ, ", "),
    ))
    nothing
end

"""
    System(exprs, vars, parameters = Variable[])

Create a system from the given `exprs`. `vars` are the given variables and determines
the variable ordering.

## Example
```julia
julia> @var x y;
julia> H = System([x^2, y^2], [y, x]);
julia> H([2, 3], 0)
2-element Array{Int64,1}:
 4
 9
```

It is also possible to declare additional variables.
```julia
julia> @var x y t a b;
julia> H = Homotopy([x^2 + a, y^2 + b^2], [x, y], [a, b]);
julia> H([2, 3], [5, 2])
2-element Array{Int64,1}:
 9
 13
```
"""
struct System
    expressions::Vector{Expression}
    variables::Vector{Variable}
    parameters::Vector{Variable}

    function System(
        exprs::Vector{Expression},
        vars::Vector{Variable},
        params::Vector{Variable},
    )
        check_vars_params(exprs, vars, params)
        new(exprs, vars, params)
    end
end

function System(
    exprs::Vector{<:Expression};
    variables::Vector{Variable} = variables(exprs),
    parameters::Vector{Variable} = Variable[],
)
    System(convert(Vector{Expression}, exprs), variables, parameters)
end

function System(
    exprs::Vector{<:Expression},
    variables::Vector{Variable},
    parameters::Vector{Variable} = Variable[],
)
    System(convert(Vector{Expression}, exprs), variables, parameters)
end

Base.hash(S::System, u::UInt64) =
    hash(S.expressions, hash(S.variables, hash(S.parameters, u)))

function Base.show(io::IO, F::System)
    if !get(io, :compact, false)
        println(io, "System of length $(length(F.expressions))")
        print(io, " $(length(F.variables)) variables: ", join(F.variables, ", "))
        if !isempty(F.parameters)
            print(io, "\n $(length(F.parameters)) parameters: ", join(F.parameters, ", "))
        end
        print(io, "\n\n")
        for i = 1:length(F)
            print(io, " ", F.expressions[i])
            i < length(F) && print(io, "\n")
        end
    else
        print(io, "[")
        for i = 1:length(F)
            print(io, F.expressions[i])
            i < length(F) && print(io, ", ")
        end
        print(io, "]")
    end
end

evaluate(F::System, x::AbstractVector) = evaluate(F.expressions, F.variables => x)
function evaluate(F::System, x::AbstractVector, p::AbstractVector)
    evaluate(F.expressions, F.variables => x, F.parameters => p)
end
(F::System)(x::AbstractVector) = evaluate(F, x)
(F::System)(x::AbstractVector, p::AbstractVector) = evaluate(F, x, p)

function Base.:(==)(F::System, G::System)
    F.expressions == G.expressions &&
    F.variables == G.variables && F.parameters == G.parameters
end

Base.size(F::System) = (length(F.expressions), length(F.variables))
Base.size(F::System, i::Integer) = size(F)[i]
Base.length(F::System) = length(F.expressions)
nvariables(F::System) = length(F.variables)
nparameters(F::System) = length(F.parameters)

variables(F::System) = F.variables
parameters(F::System) = F.parameters

degree(F::System) = degree(F.expressions, F.variables)
Base.iterate(F::System) = iterate(F.expressions)
Base.iterate(F::System, state) = iterate(F.expressions, state)


Base.push!(F::System, f::Expression) = (push!(F.expressions, f); F)
Base.append!(F::System, f::AbstractVector{Expression}) = (append!(F.expressions, f); F)

function Base.intersect(F::System, G::System)
    exprs = [F.expressions; G.expressions]
    vars = [F.variables; setdiff(G.variables, F.variables)]
    params = [F.parameters; setdiff(G.parameters, F.parameters)]
    System(exprs, vars, params)
end
function Base.intersect(F::System, G::AbstractVector{<:Expression})
    exprs = [F.expressions; G]
    vars = [F.variables; setdiff(variables(G), F.variables)]
    params = F.parameters
    System(exprs, vars, params)
end
Base.intersect(F::AbstractVector{<:Expression}, G::System) = intersect(G, F)

Base.copy(F::System) = Base.deepcopy(F)



##############
## Homotopy ##
##############
"""
    Homotopy(exprs, vars, t, parameters = Variable[])

Create a homotopy from the given `exprs`. `vars` are the given variables and determines
the variable ordering, `t` is the dedicated variable along which is "homotopied".

## Example
```julia
julia> @var x y t;
julia> H = Homotopy([x + t, y + 2t], [y, x], t);
julia> H([2, 3], 0)
2-element Array{Int64,1}:
 3
 2


julia> H([2, 3], 1)
2-element Array{Int64,1}:
 4
 4
```

It is also possible to declare additional variables.
```julia
julia> @var x y t a b;
julia> H = Homotopy([x^2 + t*a, y^2 + t*b], [x, y], t, [a, b]);
julia> H([2, 3], 1, [5, 2])
2-element Array{Int64,1}:
 9
 11
```
"""
struct Homotopy
    expressions::Vector{Expression}
    variables::Vector{Variable}
    t::Variable
    parameters::Vector{Variable}

    function Homotopy(
        exprs::Vector{Expression},
        vars::Vector{Variable},
        t::Variable,
        params::Vector{Variable},
    )
        check_vars_params(exprs, [vars; t], params)
        new(exprs, vars, t, params)
    end
end

function Homotopy(
    exprs::Vector{<:Expression},
    variables::Vector{Variable},
    t::Variable,
    parameters::Vector{Variable} = Variable[],
)
    Homotopy(convert(Vector{Expression}, exprs), variables, t, parameters)
end

function Base.show(io::IO, H::Homotopy)
    if !get(io, :compact, false)
        println(io, "Homotopy in ", H.t, " of length ", length(H.expressions))
        print(io, " $(length(H.variables)) variables: ", join(H.variables, ", "))
        if !isempty(H.parameters)
            print(io, "\n $(length(H.parameters)) parameters: ", join(H.parameters, ", "))
        end
        print(io, "\n\n")
        for i = 1:length(H)
            print(io, " ", H.expressions[i])
            i < length(H) && print(io, "\n")
        end
    else
        print(io, "[")
        for i = 1:length(H)
            print(io, H.expressions[i])
            i < length(H) && print(io, ", ")
        end
        print(io, "]")
    end
end

evaluate(H::Homotopy, x::AbstractVector, t) =
    evaluate(H.expressions, H.variables => x, H.t => t)
function evaluate(H::Homotopy, x::AbstractVector, t, p::AbstractVector)
    evaluate(H.expressions, H.variables => x, H.t => t, H.parameters => p)
end
(H::Homotopy)(x::AbstractVector, t) = evaluate(H, x, t)
(H::Homotopy)(x::AbstractVector, t, p::AbstractVector) = evaluate(H, x, t, p)

function Base.:(==)(H::Homotopy, G::Homotopy)
    H.expressions == G.expressions &&
    H.variables == G.variables && H.parameters == G.parameters
end

Base.size(H::Homotopy) = (length(H.expressions), length(H.variables))
Base.size(H::Homotopy, i::Integer) = size(H)[i]
Base.length(H::Homotopy) = length(H.expressions)
