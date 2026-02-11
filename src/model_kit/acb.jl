acb_op_call(op) = Symbol(:acb_, op_call(op), :!)

# # Arity 0
# OP_STOP
acb_op_stop!(t, m) = nothing

# # Arity 1
# OP_CB # a ^ 3
Base.@propagate_inbounds function acb_op_cb!(t, x, m)
    Arblib.sqr!(m[1], x)
    return Arblib.mul!(t, m[1], x)
end
# OP_COS # cos(a)
Base.@propagate_inbounds function acb_op_cos!(t, x, m)
    return Arblib.cos!(t, x)
end
# OP_INV # 1 / a
Base.@propagate_inbounds function acb_op_inv!(t, x, m)
    return Arblib.inv!(t, x)
end
# OP_INV_NOT_ZERO # a ≂̸ 0 ? 1 / a : a
Base.@propagate_inbounds function acb_op_inv_not_zero!(t, x, m)
    return iszero(x) ? t[] = x : Arblib.inv!(t, x)
end
# OP_INVSQR # 1 / a^2
Base.@propagate_inbounds function acb_op_invsqr!(t, x, m)
    Arblib.sqr!(m[1], x)
    return Arblib.inv!(t, m[1])
end
# OP_NEG # -a
Base.@propagate_inbounds function acb_op_neg!(t, x, m)
    return Arblib.neg!(t, x)
end
# OP_SIN # sin(a)
Base.@propagate_inbounds function acb_op_sin!(t, x, m)
    return Arblib.sin!(t, x)
end
# OP_SQR # a ^ 2
Base.@propagate_inbounds function acb_op_sqr!(t, x, m)
    return Arblib.sqr!(t, x)
end
# OP_SQRT # √(a) TODO
Base.@propagate_inbounds function acb_op_sqrt!(t, x, m)
    return Arblib.sqrt!(t, x)
end
# OP_IDENTITY # a
Base.@propagate_inbounds function acb_op_identity!(t, x, m)
    return t[] = x
end

# # Arity 2
# OP_ADD # a + b
Base.@propagate_inbounds acb_op_add!(t, x, y, m) = Arblib.add!(t, x, y)
# OP_DIV # a / b
Base.@propagate_inbounds acb_op_div!(t, x, y, m) = Arblib.div!(t, x, y)
# OP_MUL # a * b
Base.@propagate_inbounds acb_op_mul!(t, x, y, m) = Arblib.mul!(t, x, y)
# OP_SUB # a - b
Base.@propagate_inbounds acb_op_sub!(t, x, y, m) = Arblib.sub!(t, x, y)

# OP_POW_INT # a ^ p where p isa Integer
Base.@propagate_inbounds acb_op_pow_int!(t, x, k, m) = Arblib.pow!(t, x, k)
# # OP_POW # a ^ b where b isa Number

# # Arity 3
# OP_ADD3 # a + b + c
Base.@propagate_inbounds function acb_op_add3!(t, x, y, z, m)
    Arblib.add!(m[1], x, y)
    return Arblib.add!(t, m[1], z)
end
# OP_MUL3 # a * b * c TODO
Base.@propagate_inbounds function acb_op_mul3!(t, x, y, z, m)
    Arblib.mul!(m[1], x, y)
    return Arblib.mul!(t, m[1], z)
end
# OP_MULADD # a * b + c
Base.@propagate_inbounds function acb_op_muladd!(t, x, y, z, m)
    Arblib.mul!(m[1], x, y)
    return Arblib.add!(t, m[1], z)
end
# OP_MULSUB # a * b - c
Base.@propagate_inbounds function acb_op_mulsub!(t, x, y, z, m)
    Arblib.mul!(m[1], x, y)
    return Arblib.sub!(t, m[1], z)
end
# OP_SUBMUL # c - a * b
Base.@propagate_inbounds function acb_op_submul!(t, x, y, z, m)
    Arblib.mul!(m[1], x, y)
    return Arblib.sub!(t, z, m[1])
end

# # Arity 4
# OP_ADD4 # a + b + c + d
Base.@propagate_inbounds function acb_op_add4!(t, x, y, z, w, m)
    Arblib.add!(m[1], x, y)
    Arblib.add!(m[2], z, w)
    return Arblib.add!(t, m[1], m[2])
end
# OP_MUL4 # a * b * c * d TODO
Base.@propagate_inbounds function acb_op_mul4!(t, x, y, z, w, m)
    Arblib.mul!(m[1], x, y)
    Arblib.mul!(m[2], z, w)
    return Arblib.mul!(t, m[1], m[2])
end
# OP_MULMULADD # a * b + c * d
Base.@propagate_inbounds function acb_op_mulmuladd!(t, x, y, z, w, m)
    Arblib.mul!(m[1], x, y)
    Arblib.mul!(m[2], z, w)
    return Arblib.add!(t, m[1], m[2])
end
# OP_MULMULSUB # a * b - c * d
Base.@propagate_inbounds function acb_op_mulmulsub!(t, x, y, z, w, m)
    Arblib.mul!(m[1], x, y)
    Arblib.mul!(m[2], z, w)
    return Arblib.sub!(t, m[1], m[2])
end
