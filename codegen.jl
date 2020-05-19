
using Pkg
Pkg.activate(@__DIR__)

using HomotopyContinuation2

function _div_taylor_impl(N, M)
    list = ModelKit.InstructionList()
    id = push!(list, (:/, :a, :b))
    diff_map = ModelKit.DiffMap()
    for i = 1:N
        diff_map[:a, i] = :($(Symbol(:a, i)))
    end
    for i = 1:M
        diff_map[:b, i] = :($(Symbol(:b, i)))
    end

    k = max(N, M)
    dlist = ModelKit.univariate_diff!(list, k, diff_map)
    ids = [id]
    for i = 1:(k)
        push!(ids, diff_map[id, i])
    end

    quote
        $(Expr(:tuple, :a, (Symbol(:a, k) for k = 1:N)...)) = a
        $(Expr(:tuple, :b, (Symbol(:b, k) for k = 1:N)...)) = b
        $(ModelKit.to_expr(dlist))
        $(Expr(:tuple, ids...))
    end
end


@generated function div_taylor(a::NTuple{N}, b::NTuple{M}) where {N,M}
    _div_taylor_impl(N - 1, M - 1)
end

function test_div(a, b)
    da = (a, 2a, 3a, 5a)
    db = (b, -3a * b, 2.3 * b^2, 2a)
    div_taylor(da, db)
end


using BenchmarkTools

a = Ref(2.2132 + 3im)
b = Ref(1.1312 + 0.1im)

@benchmark test_div($a[], $b[])

AS = (
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
)
BS = (
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
)
us = (
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
    randn(ComplexF64, 7),
)


@benchmark ModelKit.div_taylor_vec!($us, $AS, $BS)




ModelKit._div_taylor_vec_impl(2, 2)

diff_map = ModelKit.DiffMap()
k_a = 2
for i = 1:k_a
    diff_map[:a, i] = :($(Symbol(:a, i)))
end
k_b = 2
for i = 1:k_b
    diff_map[:b, i] = :($(Symbol(:b, i)))
end

ModelKit.univariate_diff!(list, 2, diff_map)


slp = to_expr(dlist, var_map, assignements)

diff_map
