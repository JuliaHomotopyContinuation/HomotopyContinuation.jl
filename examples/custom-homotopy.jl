using HomotopyContinuation
using LinearAlgebra

struct RandomUnitaryPath{Start<:AbstractSystem,Target<:AbstractSystem} <: HomotopyContinuation.AbstractHomotopy
    straightline::StraightLineHomotopy{Start, Target}
    U::Matrix{ComplexF64}
end

function RandomUnitaryPath(start::AbstractSystem, target::AbstractSystem)
    m, n = size(start)
    # construct a random unitary matrix
    U = qr(randn(n,n) + im * randn(n,n))[1]
    RandomUnitaryPath(HomotopyContinuation.StraightLineHomotopy(start, target), U)
end

# We have to define the size
Base.size(H::RandomUnitaryPath) = size(H.straightline)

# Cache for efficient evaluation
struct RandomUnitaryPathCache{C, T1, T2} <: HomotopyContinuation.AbstractHomotopyCache
    straightline::C
    U_t::Matrix{ComplexF64}
    y::Vector{T1}
    # More temporary storage necessary to avoid allocations
    jac::Matrix{T2} # holds a jacobian
    dt::Vector{T2} # holds a derivative w.r.t. t
    U::Matrix{ComplexF64} # holds something like U
end

function HomotopyContinuation.cache(H::RandomUnitaryPath, x, t)
    U_t = copy(H.U)
    y = U_t * x
    straightline = HomotopyContinuation.cache(H.straightline, y, t)

    jac = HomotopyContinuation.jacobian(H.straightline, y, t, straightline)
    dt = jac * y
    U = copy(U_t)
    RandomUnitaryPathCache(straightline, U_t, y, jac, dt, U)
end

# Evaluation and Differentation
"""
    Ut_mul_x!(cache, U, t)

Evaluate ``U(t)x`` and return the result.
"""
function Ut_mul_x!(cache, U, x, t)
    # We start with U * (the 2x2 sin-cos block + I)
    cache.U .= U
    s, c = sin(2π*t), cos(2π*t)
    for i=1:size(U, 1)
        cache.U[i, 1] = U[i,2] * s + U[i,1] * c
        cache.U[i, 2] = U[i,2] * c - U[i,1] * s
    end
    # U(t) = cache.U * U'
    # y = U(t) * x
    mul!(cache.y, mul!(cache.U_t, cache.U, adjoint(U)), x)
end

"""
    U_dot_t_mul_x!(cache, U, t)

Evaluate ``U'(t)x`` and return the result.
"""
function U_dot_t_mul_x!(cache, U, x, t)
    # We start with U * (the derivative of the 2x2 sin-cos block + 0)
    cache.U .= zero(eltype(U))
    s, c = 2π*sin(2π*t), 2π*cos(2π*t)
    for i=1:size(U, 1)
        cache.U[i, 1] =  U[i,2] * c - U[i,1] * s
        cache.U[i, 2] = -U[i,2] * s - U[i,1] * c
    end
    # U'(t) = cache.U * U'
    # y' = U'(t) * x
    mul!(cache.y, mul!(cache.U_t, cache.U, adjoint(U)), x)
end

function HomotopyContinuation.evaluate!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    HomotopyContinuation.evaluate!(out, H.straightline, y, t, cache.straightline)
end

function HomotopyContinuation.jacobian!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    HomotopyContinuation.jacobian!(cache.jac, H.straightline, y, t, cache.straightline)
    mul!(out, cache.jac, cache.U_t) # out = J_H(y, t) * U(t)
end

function HomotopyContinuation.dt!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    # chain rule
    HomotopyContinuation.jacobian_and_dt!(cache.jac, out, H.straightline, y, t, cache.straightline)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t) # y_dot = U'(t)x
    mul!(cache.dt, cache.jac, y_dot) # dt = J_H(y, t) * y_dot
    out .+= cache.dt
end

function HomotopyContinuation.evaluate_and_jacobian!(val, jac, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    HomotopyContinuation.evaluate_and_jacobian!(val, cache.jac, H.straightline, y, t, cache.straightline)
    mul!(jac, cache.jac, cache.U_t)
end

function HomotopyContinuation.jacobian_and_dt!(jac, dt, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    HomotopyContinuation.jacobian_and_dt!(cache.jac, dt, H.straightline, y, t, cache.straightline)
    mul!(jac, cache.jac, cache.U_t) # jac = J_H(y, t) * U(t)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t) # y_dot = U'(t)x
    mul!(cache.dt, cache.jac, y_dot) # dt = J_H(y, t) * y_dot
    dt .+= cache.dt
end

# Test the implementation
# =======================
#
# We provide a test suite to test the homotopy interface. Let's use this
# to check our implementation.
@polyvar x y z

F = SPSystem([x^2*y-3x*z, z^2*x+3y^2])
G = SPSystem([z*x^2-3x*y^2, z^3*x-2x*y*z^2])
InterfaceTest.homotopy(RandomUnitaryPath(G, F))


# And now we can solve systems with our custom homotopy
solve([x^2 - y, y^3*x-x], homotopy=RandomUnitaryPath)
