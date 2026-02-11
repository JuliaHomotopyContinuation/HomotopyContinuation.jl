# More op ideas:
# eucl_sqr = (a - b) ^ 2
# inv_sqr = a ^ - 2

@enum OpType begin
    # Arity 0
    OP_STOP

    # Arity 1
    OP_CB # a ^ 3 TODO
    OP_COS # cos(a)
    OP_INV # 1 / a
    OP_INV_NOT_ZERO # a ≂̸ 0 ? 1 / a : a
    OP_INVSQR # 1 / a^2
    OP_NEG # -a
    OP_SIN # sin(a)
    OP_SQR # a ^ 2
    OP_SQRT # √(a) TODO
    OP_IDENTITY # a

    # Arity 2
    OP_ADD # a + b
    OP_DIV # a / b
    OP_MUL # a * b
    OP_SUB # a - b

    OP_POW_INT # a ^ p where p isa Integer
    # OP_POW # a ^ b where b isa Number

    # Arity 3
    OP_ADD3 # a + b + c
    OP_MUL3 # a * b * c TODO
    OP_MULADD # a * b + c
    OP_MULSUB # a * b - c
    OP_SUBMUL # c - a * b

    # Arity 4
    OP_ADD4 # a + b + c + d
    OP_MUL4 # a * b * c * d TODO
    OP_MULMULADD # a * b + c * d
    OP_MULMULSUB # a * b - c * d
end

function arity(op_type::OpType)
    return if op_type === OP_STOP
        0
    elseif op_type == OP_CB ||
            op_type == OP_COS ||
            op_type == OP_IDENTITY ||
            op_type == OP_INV ||
            op_type == OP_INV_NOT_ZERO ||
            op_type == OP_INVSQR ||
            op_type == OP_NEG ||
            op_type == OP_SIN ||
            op_type == OP_SQR ||
            op_type == OP_SQRT
        1
    elseif op_type == OP_ADD ||
            op_type == OP_DIV ||
            op_type == OP_MUL ||
            op_type == OP_SUB ||
            op_type == OP_POW_INT
        # || op_type == OP_POW
        2
    elseif op_type == OP_ADD3 ||
            op_type == OP_MUL3 ||
            op_type == OP_MULADD ||
            op_type == OP_MULSUB ||
            op_type == OP_SUBMUL
        3
    elseif op_type == OP_ADD4 ||
            op_type == OP_MUL4 ||
            op_type == OP_MULMULADD ||
            op_type == OP_MULMULSUB
        4
    else
        error("Unexpected OpType $(op_type)")
    end
end


function op_call(op_type::OpType)
    return if op_type == OP_STOP
        :op_stop

        # Arity 1
    elseif op_type == OP_CB
        :op_cb
    elseif op_type == OP_COS
        :op_cos
    elseif op_type == OP_IDENTITY
        :op_identity
    elseif op_type == OP_INV
        :op_inv
    elseif op_type == OP_INV_NOT_ZERO
        :op_inv_not_zero
    elseif op_type == OP_INVSQR
        :op_invsqr
    elseif op_type == OP_NEG
        :op_neg
    elseif op_type == OP_SIN
        :op_sin
    elseif op_type == OP_SQR
        :op_sqr
    elseif op_type == OP_SQRT
        :op_sqrt


        # Arity 2
    elseif op_type == OP_ADD
        :op_add
    elseif op_type == OP_DIV
        :op_div
    elseif op_type == OP_MUL
        :op_mul
    elseif op_type == OP_SUB
        :op_sub

    elseif op_type == OP_POW_INT
        :op_pow_int
        # elseif op_type == OP_POW
        #     :op_pow

        # Arity 3
    elseif op_type == OP_ADD3
        :op_add3
    elseif op_type == OP_MUL3
        :op_mul3
    elseif op_type == OP_MULADD
        :op_muladd
    elseif op_type == OP_MULSUB
        :op_mulsub
    elseif op_type == OP_SUBMUL
        :op_submul

        # Arity 4
    elseif op_type == OP_ADD4
        :op_add4
    elseif op_type == OP_MUL4
        :op_mul4
    elseif op_type == OP_MULMULADD
        :op_mulmuladd
    elseif op_type == OP_MULMULSUB
        :op_mulmulsub
    else
        error("Unexpected OpType $(op_type)")
    end
end

function should_use_index_not_reference(op::OpType, index::Integer)
    return (op == OP_POW_INT && index == 2)
end


# arity 0
op_stop() = nothing
# arity 1
op_cb(x) = x * x * x
@inline function op_sqr(z::Complex)
    x, y = reim(z)
    return Complex((x + y) * (x - y), (x + x) * y)
end
op_cos(x) = cos(x)
op_identity(x) = identity(x)
op_inv(x) = inv(x)
op_inv(x::Complex) = Base.FastMath.inv_fast(x)
"""
    op_inv_not_zero(x)

Invert x unless it is 0, then return 0.
"""
op_inv_not_zero(x) = ifelse(is_zero(x), x, op_inv(x))
op_invsqr(x) = op_sqr(op_inv(x))
op_neg(x) = -x
op_sin(x) = sin(x)
op_sqr(x) = x * x
op_sqrt(x) = sqrt(x)

# arity 2
op_add(a, b) = a + b
op_div(a, b) = Base.FastMath.div_fast(a, b)
op_mul(a, b) = a * b
op_sub(a, b) = a - b

op_pow_int(x, p::Integer) =
    p > 0 ? Base.power_by_squaring(x, p) : op_inv(Base.power_by_squaring(x, -p))
op_pow(x, y) = x^y

# arity 3
op_add3(x, y, z) = x + y + z
op_mul3(x, y, z) = x * y * z
op_muladd(x, y, z) = x * y + z
op_mulsub(x, y, z) = x * y - z
op_submul(x, y, z) = z - x * y

# arity 4
op_add4(a, b, c, d) = a + b + c + d
op_mul4(a, b, c, d) = a * b * c * d
op_mulmuladd(a, b, c, d) = a * b + c * d
op_mulmulsub(a, b, c, d) = a * b - c * d
