########################
## Lift to Type Level ##
########################

abstract type TExpression end
struct TOperation{F,A1,A2} <: TExpression end

TExpression(x) = convert(TExpression, x)
Base.convert(::Type{TExpression}, v::Variable) = v.name
Base.convert(::Type{TExpression}, x::Constant) = x.value
function Base.convert(::Type{TExpression}, op::Operation)
    if length(op.args) == 1
        TOperation{op.func,TExpression(op.args[1]),Nothing}()
    else
        TOperation{op.func,TExpression(op.args[1]),TExpression(op.args[2])}()
    end
end

Base.convert(::Type{Expression}, ::TOperation{F,A1,Nothing}) where {F,A1} =
    Operation(F, Expression[Expression(A1)])
Base.convert(::Type{Expression}, ::TOperation{F,A1,A2}) where {F,A1,A2} =
    Operation(F, Expression[Expression(A1), Expression(A2)])

#############
## TSystem ##
#############
struct TSystem{TE,V,P} end

type_level(sys::System) = typeof(TSystem(sys.expressions, sys.variables, sys.parameters))
function TSystem(
    exprs::Vector{<:Expression},
    var_order::AbstractVector{<:Variable},
    param_order::AbstractVector{<:Variable} = Variable[],
)
    V = tuple((var.name for var in var_order)...)
    if !isempty(param_order)
        P = tuple((var.name for var in param_order)...)
    else
        P = Nothing
    end
    TE = tuple(convert.(TExpression, exprs)...)
    TSystem{TE,V,P}()
end

Base.show(io::IO, ::Type{T}) where {T<:TSystem} = show_info(io, T)
function Base.show(io::IO, TS::TSystem)
    show_info(io, typeof(TS))
    print(io, "()")
end
function show_info(io::IO, ::Type{TSystem{TE,V,P}}) where {TE,V,P}
    n = length(TE)
    m = length(V)
    mp = (P == Nothing ? 0 : length(P))
    print(io, "TSystem{$n,$m,$mp,#$(hash(TE))}")
end


interpret(TS::TSystem) = interpret(typeof(TS))
function interpret(::Type{TSystem{TE,V,P}}) where {TE,V,P}
    exprs = [Expression(e) for e in TE]
    vars = [Variable(e) for e in V]
    if P == Nothing
        params = Variable[]
    else
        params = [Variable(v) for v in P]
    end
    System(exprs, vars, params)
end

###############
## THomotopy ##
###############

struct THomotopy{TE,V,T,P} end

type_level(sys::Homotopy) =
    typeof(THomotopy(sys.expressions, sys.variables, sys.t, sys.parameters))
function THomotopy(
    exprs::Vector{<:Expression},
    var_order::AbstractVector{<:Variable},
    t::Variable,
    param_order::AbstractVector{<:Variable} = Variable[],
)
    TE = tuple(convert.(TExpression, exprs)...)
    V = tuple((var.name for var in var_order)...)
    T = t.name
    P = isempty(param_order) ? Nothing : tuple((var.name for var in param_order)...)
    THomotopy{TE,V,T,P}()
end

Base.show(io::IO, ::Type{T}) where {T<:THomotopy} = show_info(io, T)
function Base.show(io::IO, TS::THomotopy)
    show_info(io, typeof(TS))
    print(io, "()")
end
function show_info(io::IO, ::Type{THomotopy{TE,V,T,P}}) where {TE,V,T,P}
    n = length(TE)
    m = length(V)
    mp = (P == Nothing ? 0 : length(P))
    print(io, "THomotopy{$n,$m,$mp,#$(hash(TE))}")
end

interpret(TS::THomotopy) = interpret(typeof(TS))
function interpret(::Type{THomotopy{TE,V,T,P}}) where {TE,V,T,P}
    exprs = [Expression(e) for e in TE]
    vars = [Variable(e) for e in V]
    t = Variable(T)
    params = P == Nothing ? Variable[] : [Variable(v) for v in P]
    Homotopy(exprs, vars, t, params)
end

###########################
## CODEGEN + COMPILATION ##
###########################

struct InstructionId
    name::Symbol
end
name(id::InstructionId) = id.name

Base.:(==)(x::InstructionId, y::InstructionId) = x.name == y.name
Base.show(io::IO, v::InstructionId) = print(io, v.name)
Base.convert(::Type{Expr}, x::InstructionId) = x.name

const SimpleExpression = Union{InstructionId,Constant,Variable}

struct Instruction
    func::Union{Nothing,Symbol} # module, function name
    args::Vector{SimpleExpression}
end

Instruction(f::Union{Nothing,Symbol}, a::SimpleExpression) = Instruction(f, [a])
Instruction(x::SimpleExpression) = Instruction(nothing, x)

Base.hash(I::Instruction, h::UInt) = foldr(hash, I.args, init = hash(I.func, h))
Base.:(==)(a::Instruction, b::Instruction) = a.func == b.func && a.args == b.args

function Base.convert(::Type{Expr}, op::Instruction)
    if op.func === nothing && length(op.args) == 1
        convert(Expr, op.args[1])
    else
        Expr(:call, op.func, convert.(Expr, op.args)...)
    end
end
Base.show(io::IO, op::Instruction) = print(io, convert(Expr, op))


struct InstructionList
    instructions::OrderedDict{Instruction,InstructionId}
    var::Symbol
    n::Base.RefValue{Int}
end

function InstructionList(; var::Symbol = :ι, n::Base.RefValue{Int} = Ref(0))
    instructions = OrderedDict{Instruction,InstructionId}()
    InstructionList(instructions, var, n)
end

function Base.push!(v::InstructionList, i::Instruction; id_var::Symbol = v.var)
    if haskey(v.instructions, i)
        v.instructions[i]
    else
        id = InstructionId(Symbol(id_var, v.n[] += 1))
        push!(v.instructions, i => id)
        id
    end
end

Base.push!(list::InstructionList, x::SimpleExpression) = push!(list, Instruction(x))

function Base.push!(list::InstructionList, op::Operation)
    if op.func == :^ &&
       op.args[2] isa Constant && op.args[2].value isa Integer && op.args[2].value > 1
        n = op.args[2].value
        if n == 2
            i = _push!(list, op.args[1])
            push!(list, Instruction(:*, [i, i]))
        else
            n₁ = div(n, 2) + 1 # make sure n₁ > n₂
            n₂ = n - n₁
            id1 = push!(list, Operation(:^, op.args[1], Constant(n₁)))
            if n₂ == 1
                push!(list, Instruction(:*, [id1, _push!(list, op.args[1])]))
            else
                id2 = push!(list, Operation(:^, op.args[1], Constant(n₂)))
                push!(list, Instruction(:*, [id1, id2]))
            end
        end
    else
        push!(list, Instruction(op.func, map(arg -> _push!(list, arg), op.args)))
    end
end
_push!(list::InstructionList, op::Operation) = push!(list, op)
_push!(list::InstructionList, op::Expression) = op

Base.length(v::InstructionList) = length(v.instructions)
Base.iterate(v::InstructionList) = iterate(v.instructions)
Base.iterate(v::InstructionList, state) = iterate(v.instructions, state)
Base.eltype(v::Type{InstructionList}) = eltype(v.instructions)

function Base.show(io::IO, ::MIME"text/plain", list::InstructionList)
    println(io, "InstructionList:")
    for (instr, id) in list
        println(io, convert(Expr, id), " = ", convert(Expr, instr))
    end
end

function Base.convert(::Type{Expr}, list::InstructionList)
    Expr(:block, map(list) do (instr, id)
        :($(convert(Expr, id)) = $(convert(Expr, instr)))
    end...)
end


## CODEGEN
struct Compiled{T1<:Union{TSystem,THomotopy},T2<:Union{System,Homotopy}}
    obj::T2
end
Compiled(obj::T) where {T<:Union{System,Homotopy}} = Compiled{type_level(obj),T}(obj)

Base.size(C::Compiled) = size(C.obj)
Base.size(C::Compiled, i::Integer) = size(C.obj, i)
Base.length(C::Compiled) = length(C.obj)

const CompiledSystem{T<:TSystem} = Compiled{T,System}
const CompiledHomotopy{T<:THomotopy} = Compiled{T,Homotopy}

compile(F::System) = Compiled(F)
compile(H::Homotopy) = Compiled(H)

interpreted(C::Compiled) = C.obj

function Base.show(io::IO, C::Compiled{T}) where {T}
    println(io, "Compiled{", T, "}:")
    show(io, C.obj)
end

#################
## evaluations ##
#################
const TExpr = Union{TSystem,THomotopy}

make_indexing(F::System) = make_indexing(F.variables, nothing, F.parameters)
make_indexing(H::Homotopy) = make_indexing(H.variables, H.t, H.parameters)
function make_indexing(vars, t, params)
    lhs = Any[convert(Expr, v) for v in vars]
    rhs = Any[:(x[$i]) for i = 1:length(vars)]
    if t !== nothing
        push!(lhs, convert(Expr, t))
        push!(rhs, :t)
    end
    if params !== nothing
        append!(lhs, convert.(Expr, params))
        for i = 1:length(params)
            push!(rhs, :(p[$i]))
        end
    end
    :($(Expr(:tuple, lhs...)) = $(Expr(:tuple, rhs...)))
end

interpret(T::Type{<:TSystem}) = System(T)
interpret(T::Type{<:THomotopy}) = Homotopy(T)

function jacobian!(U::AbstractMatrix, S::Compiled, args...)
    evaluate_and_jacobian!(nothing, U, S, args...)
    U
end
function evaluate_and_jacobian(S::Compiled, args...)
    evaluate(S, args...), jacobian(S, args...)
end

function _evaluate_impl(::Type{T}) where {T<:TExpr}
    I = interpret(T)
    quote
        let $(make_indexing(I))
            $(begin
                list = InstructionList()
                ids = [convert(Expr, push!(list, fi)) for fi in I.expressions]
                quote
                    $(convert(Expr, list))
                    @SVector $(Expr(:vect, ids...))
                end
            end)
        end
    end
end

@generated function evaluate(F::CompiledSystem{T}, x::AbstractVector, p = nothing) where {T}
    _evaluate_impl(T)
end
(F::CompiledSystem)(x, p = nothing) = evaluate(F, x, p)
@generated function evaluate(
    F::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p = nothing,
) where {T}
    _evaluate_impl(T)
end
(H::CompiledHomotopy)(x, t, p = nothing) = evaluate(H, x, t, p)

function _evaluate!_impl(::Type{T}) where {T<:TExpr}
    u = gensym(:u)
    I = interpret(T)
    quote
        let $u = u
            let $(make_indexing(I))
                $(begin
                    list = InstructionList()
                    ids = [convert(Expr, push!(list, fi)) for fi in I.expressions]
                    quote
                        $(convert(Expr, list))
                        $(map(i -> :($u[$i] = $(ids[i])), 1:length(ids))...)
                    end
                end)
            end
        end
        u
    end
end

@generated function evaluate!(
    u::AbstractVector,
    F::CompiledSystem{T},
    x::AbstractVector,
    p = nothing,
) where {T}
    _evaluate!_impl(T)
end
@generated function evaluate!(
    u::AbstractVector,
    F::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p = nothing,
) where {T}
    _evaluate!_impl(T)
end

function _evaluate_and_jacobian!_impl(::Type{T}) where {T<:TExpr}
    u, U = gensym(:u), gensym(:U)
    I = interpret(T)
    n, m = size(I)
    quote
        let ($u, $U) = (u, U)
            let $(make_indexing(I))
                $(begin
                    list = InstructionList()
                    eval_ids = Vector{Any}(undef, n)
                    jac_ids = Matrix{Any}(undef, n, m)
                    for i = 1:n
                        fᵢ = I.expressions[i]
                        eval_ids[i] = convert(Expr, push!(list, fᵢ))
                        for j = 1:m
                            jac_ids[i, j] = convert(Expr, push!(list, I.jacobian[i, j]))
                        end
                    end
                    quote
                        $(convert(Expr, list))
                        if !($u isa Nothing)
                            $(map(i -> :($u[$i] = $(eval_ids[i])), 1:n)...)
                        end
                        $(vec([:($U[$i, $j] = $(jac_ids[i, j])) for i = 1:n, j = 1:m])...)
                        nothing
                    end
                end)
            end
        end
    end
end

@generated function evaluate_and_jacobian!(
    u::Union{Nothing,AbstractVector},
    U::AbstractMatrix,
    F::CompiledSystem{T},
    x::AbstractVector,
    p::Union{Nothing,AbstractVector} = nothing,
) where {T}
    _evaluate_and_jacobian!_impl(T)
end
@generated function evaluate_and_jacobian!(
    u::Union{Nothing,AbstractVector},
    U::AbstractMatrix,
    F::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p::Union{Nothing,AbstractVector} = nothing,
) where {T}
    _evaluate_and_jacobian!_impl(T)
end

function _jacobian_impl(::Type{T}) where {T<:TExpr}
    I = interpret(T)
    n, m = size(I)
    quote
        let $(make_indexing(I))
            $(begin
                list = InstructionList()
                jac_ids = map(e -> convert(Expr, push!(list, e)), I.jacobian)
                quote
                    $(convert(Expr, list))
                    @SMatrix $(Expr(:vcat, (Expr(:row, jac_ids[i, :]...) for i = 1:n)...))
                end
            end)
        end
    end
end

@generated function jacobian(F::CompiledSystem{T}, x, p = nothing) where {T}
    _jacobian_impl(T)
end

@generated function jacobian(F::CompiledHomotopy{T}, x, t, p = nothing) where {T}
    _jacobian_impl(T)
end

function _dt_impl(::Type{T}) where {T<:THomotopy}
    I = interpret(T)
    quote
        let $(make_indexing(I))
            $(begin
                list = InstructionList()
                ids = map(e -> convert(Expr, push!(list, e)), I.dt)
                quote
                    $(convert(Expr, list))
                    @SVector $(Expr(:vect, ids...))
                end
            end)
        end
    end
end

@generated function dt(H::CompiledHomotopy{T}, x::AbstractVector, t, p = nothing) where {T}
    _dt_impl(T)
end

function _dt!_impl(::Type{T}) where {T<:THomotopy}
    u = gensym(:u)
    I = interpret(T)
    quote
        let $u = u
            let $(make_indexing(I))
                $(begin
                    list = InstructionList()
                    ids = map(e -> convert(Expr, push!(list, e)), I.dt)
                    quote
                        $(convert(Expr, list))
                        $((:($u[$i] = $(ids[i])) for i = 1:length(ids))...)
                    end
                end)
            end
        end
        u
    end
end

@generated function dt!(
    u::AbstractVector,
    H::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p = nothing,
) where {T}
    _dt!_impl(T)
end

####

function _jacobian_and_dt!_impl(::Type{T}) where {T<:THomotopy}
    u, U = gensym(:u), gensym(:U)
    I = interpret(T)
    n, m = size(I)
    quote
        let ($u, $U) = (u, U)
            let $(make_indexing(I))
                $(begin
                    list = InstructionList()
                    dt_ids = Vector{Any}(undef, n)
                    jac_ids = Matrix{Any}(undef, n, m)
                    for i = 1:n
                        fᵢ = I.expressions[i]
                        dt_ids[i] = convert(Expr, push!(list, I.dt[i]))
                        for j = 1:m
                            jac_ids[i, j] = convert(Expr, push!(list, I.jacobian[i, j]))
                        end
                    end
                    quote
                        $(convert(Expr, list))
                        if !($u isa Nothing)
                            $(map(i -> :($u[$i] = $(dt_ids[i])), 1:n)...)
                        end
                        $(vec([:($U[$i, $j] = $(jac_ids[i, j])) for i = 1:n, j = 1:m])...)
                        nothing
                    end
                end)
            end
        end
    end
end

@generated function jacobian_and_dt!(
    u::AbstractVector,
    U::AbstractMatrix,
    H::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p = nothing,
) where {T}
    _jacobian_and_dt!_impl(T)
end

function jacobian_and_dt(
    H::CompiledHomotopy{T},
    x::AbstractVector,
    t,
    p = nothing,
) where {T}
    dt(H, x, t, p), jacobian(H, x, t, p)
end
