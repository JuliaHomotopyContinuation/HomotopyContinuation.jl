struct TruncatedTaylorSeries{N,T}
    val::NTuple{N,T}
end
TruncatedTaylorSeries(v::AbstractVector) = TruncatedTaylorSeries(tuple(v...))
TruncatedTaylorSeries(v::Number) = TruncatedTaylorSeries((v,))
TruncatedTaylorSeries(v::Tuple) = TruncatedTaylorSeries(promote(v...))

Base.show(io::IO, x::TruncatedTaylorSeries) = show(io, x.val)
Base.eltype(::TruncatedTaylorSeries{N,T}) where {N,T} = T
Base.length(x::TruncatedTaylorSeries) = length(x.val)
Base.iterate(x::TruncatedTaylorSeries{N,T}) where {N,T} = iterate(x.val)
Base.iterate(x::TruncatedTaylorSeries{N,T}, s) where {N,T} = iterate(x.val, s)
Base.getindex(x::TruncatedTaylorSeries{N,T}, k) where {N,T} = getindex(x.val, k + 1)
Base.zero(::TruncatedTaylorSeries{N,T}) where {N,T} =
    convert(TruncatedTaylorSeries{N,T}, zero(T))
Base.zero(::Type{TruncatedTaylorSeries{N,T}}) where {N,T} =
    convert(TruncatedTaylorSeries{N,T}, zero(T))
Base.:(==)(tx::TruncatedTaylorSeries, ty::TruncatedTaylorSeries) = tx.val == ty.val

function Base.convert(
    ::Type{TruncatedTaylorSeries{N,T}},
    x::TruncatedTaylorSeries{K1,T1},
) where {N,T,K1,T1}
    TruncatedTaylorSeries(ntuple(Val(N)) do i
        if i <= length(x.val)
            convert(T, x.val[i])
        else
            zero(T)
        end
    end)
end
function Base.convert(::Type{TruncatedTaylorSeries{N,T}}, x::T1) where {N,T,T1<:Number}
    TruncatedTaylorSeries(ntuple(Val(N)) do i
        if i == 1
            convert(T, x)
        else
            zero(T)
        end
    end)
end
function Base.convert(::Type{TruncatedTaylorSeries{N,T}}, x::Tuple) where {N,T}
    convert(TruncatedTaylorSeries{N,T}, TruncatedTaylorSeries(x))
end


Base.:*(a::Number, x::TruncatedTaylorSeries) = TruncatedTaylorSeries(a .* x.val)
Base.:*(x::TruncatedTaylorSeries, a::Number) = TruncatedTaylorSeries(a .* x.val)

for i = 1:5, j = 1:5
    i != j || continue
    @eval Base.:+(
        x::TruncatedTaylorSeries{$i,T1},
        y::TruncatedTaylorSeries{$j,T2},
    ) where {T1,T2} = TruncatedTaylorSeries(
        convert(TruncatedTaylorSeries{$(max(i, j)),T1}, x).val .+
        convert(TruncatedTaylorSeries{$(max(i, j)),T2}, y).val,
    )
    @eval Base.:-(
        x::TruncatedTaylorSeries{$i,T1},
        y::TruncatedTaylorSeries{$j,T2},
    ) where {T1,T2} = TruncatedTaylorSeries(
        convert(TruncatedTaylorSeries{$(max(i, j)),T1}, x).val .-
        convert(TruncatedTaylorSeries{$(max(i, j)),T2}, y).val,
    )
end

Base.:+(x::TruncatedTaylorSeries{N}, y::TruncatedTaylorSeries{N}) where {N} =
    TruncatedTaylorSeries(x.val .+ y.val)
Base.:-(x::TruncatedTaylorSeries{N}, y::TruncatedTaylorSeries{N}) where {N} =
    TruncatedTaylorSeries(x.val .- y.val)


expression(x, ε = Variable(:ε)) = x
function expression(x::TruncatedTaylorSeries, ε = Variable(:ε))
    ex = zero(Expression)
    for (k, xk) in enumerate(x)
        ex += xk * ε^(k - 1)
    end
    ex
end

"""
    TaylorVector{N,T}

A data structure representing a vector of Taylor series with `N` terms.
Each element is an `NTuple{N,T}`.

    TaylorVector(T, N::Integer, n::Integer)

Create a vector of `n` `NTuple{N,T}`s.
"""
struct TaylorVector{N,T} <: AbstractVector{TruncatedTaylorSeries{N,T}}
    data::Matrix{T}
end
TaylorVector{N}(data::Matrix{T}) where {N,T} = TaylorVector{N,T}(data)
TaylorVector{N}(TV::TaylorVector{M,T}) where {N,M,T} = TaylorVector{N,T}(TV.data)
TaylorVector{N}(T, n::Integer) where {N} = TaylorVector{N}(zeros(T, N, n))
TaylorVector{N,T}(::UndefInitializer, dims::Tuple{Int}) where {N,T} =
    TaylorVector{N,T}(similar(Matrix{T}, N, dims[1]))


function Base.show(io::IO, ::MIME"text/plain", X::TaylorVector)
    summary(io, X)
    isempty(X) && return
    println(io, ":")
    Base.print_array(io, map(i -> X[i], 1:length(X)))
end

Base.length(TV::TaylorVector) = size(TV.data, 2)
Base.size(TV::TaylorVector) = (length(TV),)
Base.eltype(::Type{TaylorVector{N,T}}) where {N,T} = TruncatedTaylorSeries{N,T}
Base.IndexStyle(::TaylorVector) = Base.IndexLinear()

Base.fill!(v::TaylorVector, x) = fill!(v.data, x)

"""
    vectors(TV::TaylorVec{N})

Return the Taylor series as `N` seperate vectors.
"""
@generated function vectors(TV::TaylorVector{N}) where {N}
    Expr(:tuple, (:(view(TV.data, $k, :)) for k = 1:N)...)
end

@generated function Base.getindex(TV::TaylorVector{N}, i::Integer) where {N}
    quote
        Base.@_propagate_inbounds_meta
        x = TV.data
        TruncatedTaylorSeries($(Expr(:tuple, (:(x[$k, i]) for k = 1:N)...)))
    end
end
Base.getindex(TV::TaylorVector{N,T}, i::Integer, j::Integer) where {N,T} =
    j <= N ? getindex(TV.data, j, i) : zero(T)

Base.getindex(
    TV::SubArray{<:TruncatedTaylorSeries{N,T}},
    i::Integer,
    j::Integer,
) where {N,T} = j <= N ? TV[i][j-1] : zero(T)

function Base.setindex!(TV::TaylorVector{N,T}, x, i::Integer) where {N,T}
    setindex!(TV, convert(TruncatedTaylorSeries{N,T}, x), i)
end
@generated function Base.setindex!(
    TV::TaylorVector{N,T},
    x::TruncatedTaylorSeries{N,T},
    i::Integer,
) where {N,T}
    quote
        Base.@_propagate_inbounds_meta
        d = TV.data
        $((:(d[$k, i] = x.val[$k]) for k = 1:N)...)
        x
    end
end
Base.setindex!(TV::TaylorVector, x, i::Integer, j::Integer) = setindex!(TV.data, x, j, i)

function Base.setindex!(
    TV::SubArray{<:TruncatedTaylorSeries{N,T}},
    x,
    i::Integer,
    j::Integer,
) where {N,T}

    if j <= N
        d = StaticArrays.MVector{N}(TV[i].val)
        d[j] = x

        TV[i] = StaticArrays.SVector{N}(d).data
    end
    TV
end

# HIGHER ORDER DIFFERENTIATION
# This follows Chapter 13 of Griewank & Walther - Evaluating derivatives


struct DiffMap
    D::Dict{Tuple{Any,Int},Any}
end
DiffMap() = DiffMap(Dict{Tuple{Any,Int},Any}())
#
Base.getindex(D::DiffMap, var, ∂i) = get(D.D, (var, ∂i), nothing)
Base.setindex!(D::DiffMap, val, var, ∂i) = D.D[(var, ∂i)] = val
Base.setindex!(D::DiffMap, ::Nothing, var, ∂i) = nothing


# Taylor
function untuple(varname, dx)
    :($(Expr(:tuple, (Symbol(varname, k) for k = 0:dx)...)) = $varname)
end
function taylor_tuple(ids)
    tuple = Expr(:tuple, [isnothing(id) ? :(zero(x0)) : to_expr_arg(id) for id in ids]...)
    :(TruncatedTaylorSeries($tuple))
end

function taylor_impl(f!, K::Int, dx::Int)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    ids = f!(list, D)
    quote
        Base.@_inline_meta
        $(untuple(:x, dx))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end
function taylor_impl(f!, K::Int, dx::Int, dy::Int)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    for k = 0:dy
        D[:y, k] = Symbol(:y, k)
    end
    ids = f!(list, D)
    quote
        Base.@_inline_meta
        $(untuple(:x, dx))
        $(untuple(:y, dy))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end
function taylor_impl(f!, K::Int, dx::Int, dy::Int, dz::Int)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    for k = 0:dy
        D[:y, k] = Symbol(:y, k)
    end
    for k = 0:dz
        D[:z, k] = Symbol(:z, k)
    end

    ids = f!(list, D)
    quote
        Base.@_inline_meta
        $(untuple(:x, dx))
        $(untuple(:y, dy))
        $(untuple(:z, dz))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end
function taylor_impl(f!, K::Int, dx::Int, dy::Int, dz::Int, dw::Int)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    for k = 0:dy
        D[:y, k] = Symbol(:y, k)
    end
    for k = 0:dz
        D[:z, k] = Symbol(:z, k)
    end
    for k = 0:dw
        D[:w, k] = Symbol(:w, k)
    end

    ids = f!(list, D)
    quote
        Base.@_inline_meta
        $(untuple(:x, dx))
        $(untuple(:y, dy))
        $(untuple(:z, dz))
        $(untuple(:w, dw))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end


taylor_op_call(op::OpType) = Symbol(:taylor_, op_call(op))

# Generate fallbacks to convert every input to a tuple
truncated_taylor_series(t::N) where {N<:Number} = TruncatedTaylorSeries((t,))
truncated_taylor_series(t) = TruncatedTaylorSeries(t)
truncated_taylor_series(t::T) where {T<:TruncatedTaylorSeries} = t
for op in instances(OpType)
    f = taylor_op_call(op)
    if arity(op) === 1
        @eval @inline $(f)(::Val{K}, x::X) where {K,X} =
            $f(Val(K), truncated_taylor_series(x))
    elseif arity(op) == 2
        if op == OP_POW_INT
            @eval $f(::Val{K}, x::X, r::Integer) where {K,X} =
                $f(Val(K), truncated_taylor_series(x), r)
        else
            @eval $f(::Val{K}, x::X, y::Y) where {K,X,Y} =
                $f(Val(K), truncated_taylor_series(x), truncated_taylor_series(y))
        end
    elseif arity(op) == 3
        @eval @inline $f(::Val{K}, x::X, y::Y, z::Z) where {K,X,Y,Z} = $f(
            Val(K),
            truncated_taylor_series(x),
            truncated_taylor_series(y),
            truncated_taylor_series(z),
        )
    elseif arity(op) == 4
        @eval @inline $f(::Val{K}, x::X, y::Y, z::Z, w::W) where {K,X,Y,Z,W} = $f(
            Val(K),
            truncated_taylor_series(x),
            truncated_taylor_series(y),
            truncated_taylor_series(z),
            truncated_taylor_series(w),
        )
    end
end



# OP_STOP
function taylor_op_stop(V::Val{K}) where {K}
    nothing
end

#
# OP_CB # a ^ 3 TODO
function taylor_op_cb(V::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_op_mul(V, taylor_op_sqr(V, x), x)
end
# OP_COS and OP_SIN
function taylor_sin_helper!(list, D, k)
    k == 0 && return :s₀
    s_k = nothing
    for j = 1:k
        s_k =
            muladd!(list, mul!(list, j, D[:x, j]), taylor_cos_helper!(list, D, k - j), s_k)
    end
    div!(list, s_k, k)
end
function taylor_cos_helper!(list, D, k)
    k == 0 && return :c₀
    c_k = nothing
    for j = 1:k
        c_k =
            submul!(list, mul!(list, j, D[:x, j]), taylor_sin_helper!(list, D, k - j), c_k)
    end
    div!(list, c_k, k)
end
function taylor_sin_impl(K, M)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:M-1
        D[:x, k] = Symbol(:x, k)
    end

    ids = Any[:s₀]
    for k = 1:K
        push!(ids, taylor_sin_helper!(list, D, k))
    end

    quote
        Base.@_inline_meta
        $(untuple(:x, M - 1))
        $(K > 0 ? :((s₀, c₀) = sincos(x0)) : :(s₀ = sin(x0)))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end
function taylor_cos_impl(K, M)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:M-1
        D[:x, k] = Symbol(:x, k)
    end

    ids = Any[:c₀]
    for k = 1:K
        push!(ids, taylor_cos_helper!(list, D, k))
    end

    quote
        Base.@_inline_meta
        $(untuple(:x, M - 1))
        $(K > 0 ? :((s₀, c₀) = sincos(x0)) : :(c₀ = cos(x0)))
        $(to_julia_expr(list))
        $(taylor_tuple(ids))
    end
end
@generated function taylor_op_sin(::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_sin_impl(K, M)
end
@generated function taylor_op_cos(::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_cos_impl(K, M)
end

# OP_INV # 1 / a
@generated function taylor_op_inv(V::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_impl(K, M - 1) do list, D
        ids = []
        for k = 0:K
            s = nothing
            for j = 0:(k-1)
                s = muladd!(list, ids[j+1], D[:x, k-j], s)
            end
            push!(ids, div!(list, sub!(list, k == 0 ? 1 : nothing, s), D[:x, 0]))
        end
        ids
    end
end
# OP_INV_NOT_ZERO # a ≂̸ 0 ? 1 / a : a
function taylor_op_inv_not_zero(V::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_op_div(V, 1, x)
end
# OP_INVSQR # 1 / a^2
function taylor_op_invsqr(V::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_op_div(V, 1, taylor_op_sqr(V, x))
end
# OP_NEG # -a
@generated function taylor_op_neg(::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_impl(K, M - 1) do list, D
        [neg!(list, D[:x, k]) for k = 0:K]
    end
end
# OP_SQR # a ^ 2
function taylor_op_sqr_impl(K, M)
    taylor_impl(K, M - 1) do list, D
        map(0:K) do k
            # Compute ∑_{j=0}^k u_j U_{k - j} =
            #         2(∑_{j=0}^div(k - 1,2) u_j U_{k - j}) +
            #         iseven(k) * u_{k-1}^2
            k == 0 && return sqr!(list, D[:x, 0])
            w_k = nothing
            for j = 0:div(k - 1, 2)
                w_k = muladd!(list, D[:x, j], D[:x, k-j], w_k)
            end
            if iseven(k)
                muladd!(list, 2, w_k, sqr!(list, D[:x, div(k, 2)]))
            else
                add!(list, w_k, w_k)
            end
        end
    end
end
@generated function taylor_op_sqr(::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_op_sqr_impl(K, M)
end
# OP_SQRT # √(a)
# TODO VERIFY CORRECTNESS
@generated function taylor_op_sqrt(::Val{K}, x::TruncatedTaylorSeries{M}) where {K,M}
    taylor_impl(K, M - 1) do list, D
        v₀ = add_op!(list, OP_SQRT, D[:x, 0])
        d = mul!(list, 2, v₀)
        ids = []
        push!(ids, v₀)
        if K >= 1
            v₁ = div!(list, D[:x, 1], d)
            push!(ids, v₁)
        end
        for k = 2:K
            s = nothing
            for j = 1:(k-1)
                s = muladd!(list, ids[j+1], ids[k-j+1], s)
            end
            push!(ids, div!(list, sub!(list, D[:x, k], s), d))
        end
        ids
    end
end
# OP_IDENTITY # a
function taylor_op_identity(::Val{K}, x::TruncatedTaylorSeries{M,T}) where {K,M,T}
    convert(TruncatedTaylorSeries{K + 1,T}, x)
end
#
# OP_ADD # a + b
@generated function taylor_op_add(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        [add!(list, D[:x, k], D[:y, k]) for k = 0:K]
    end
end
# OP_DIV # a / b
@generated function taylor_op_div(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        ids = []
        for k = 0:K
            s = nothing
            for j = 0:(k-1)
                s = muladd!(list, ids[j+1], D[:y, k-j], s)
            end
            push!(ids, div!(list, sub!(list, D[:x, k], s), D[:y, 0]))
        end
        ids
    end
end
# OP_MUL # a * b
@generated function taylor_op_mul(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        map(0:K) do k
            c_k = nothing
            for j = 0:k
                c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end
# OP_SUB # a - b
@generated function taylor_op_sub(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        [sub!(list, D[:x, k], D[:y, k]) for k = 0:K]
    end
end
#
# OP_POW_INT # a ^ p where p isa Integer

function taylor_op_pow_int_impl(K, dx)
    D = DiffMap()
    list = IntermediateRepresentation()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    w = Any[:w₀]
    for k = 1:K
        s = nothing
        for j = 1:k
            ũ_j = mul!(list, j, D[:x, j])
            w_kj = k == j ? :w₀ : w[k-j+1]
            s = muladd!(list, w_kj, ũ_j, s)
        end
        s = mul!(list, :r, s)

        t = nothing
        for j = 1:k-1
            u_kj = D[:x, k-j]
            w̃_j = mul!(list, j, w[j+1])
            t = muladd!(list, u_kj, w̃_j, t)
        end
        w̃_k = mul!(list, :u₀_inv, sub!(list, s, t))
        w_k = div!(list, w̃_k, k)
        push!(w, w_k)
    end
    quote
        Base.@_inline_meta
        $(untuple(:x, dx))
        iszero(x0) && return $(taylor_tuple([nothing for _ = 0:K]))
        w₀ = op_pow_int(x0, r)
        $((K > 0 ? (:(u₀_inv = op_inv(x0)),) : ())...)
        $(to_julia_expr(list))
        $(taylor_tuple(w))
    end
end
@generated function taylor_op_pow_int(
    ::Val{K},
    x::TruncatedTaylorSeries{M,T},
    r::I,
) where {K,M,T,I<:Integer}
    taylor_op_pow_int_impl(K, M - 1)
end
# OP_POW # a ^ b  TODO


# OP_ADD3 # a + b + c
function taylor_op_add3(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
) where {K,M,N,O}
    taylor_op_add(V, taylor_op_add(V, x, y), z)
end
# OP_MUL3 # a * b * c
function taylor_op_mul3(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
) where {K,M,N,O}
    taylor_op_mul(V, taylor_op_mul(V, x, y), z)
end
# OP_MULADD # a * b + c
@generated function taylor_op_muladd(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            for j = 0:k
                c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end
# OP_MULSUB # a * b - c
@generated function taylor_op_mulsub(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            neg = false
            for j = 0:k
                if neg
                    c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
                else
                    c_k = mulsub!(list, D[:x, j], D[:y, k-j], c_k)
                    neg = true
                end
            end
            neg ? c_k : neg!(list, c_k)
        end
    end
end
# OP_SUBMUL # c - a * b
@generated function taylor_op_submul(
    ::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            for j = 0:k
                c_k = submul!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end
#

# OP_ADD4 # a + b + c + d
function taylor_op_add4(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
    w::TruncatedTaylorSeries{P},
) where {K,M,N,O,P}
    taylor_op_add(V, taylor_op_add(V, x, y), taylor_op_add(V, z, w))
end
# OP_MUL4 # a * b * c * d
function taylor_op_mul4(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
    w::TruncatedTaylorSeries{P},
) where {K,M,N,O,P}
    taylor_op_mul(V, taylor_op_mul(V, x, y), taylor_op_mul(V, z, w))
end
# OP_MULMULADD # a * b + c * d
function taylor_op_mulmuladd(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
    w::TruncatedTaylorSeries{P},
) where {K,M,N,O,P}
    taylor_op_add(V, taylor_op_mul(V, x, y), taylor_op_mul(V, z, w))
end
# OP_MULMULSUB # a * b - c * d
function taylor_op_mulmulsub(
    V::Val{K},
    x::TruncatedTaylorSeries{M},
    y::TruncatedTaylorSeries{N},
    z::TruncatedTaylorSeries{O},
    w::TruncatedTaylorSeries{P},
) where {K,M,N,O,P}
    taylor_op_sub(V, taylor_op_mul(V, x, y), taylor_op_mul(V, z, w))
end
