using HomotopyContinuation

struct RandomUnitaryPath{Start<:Systems.AbstractSystem,Target<:Systems.AbstractSystem} <: Homotopies.AbstractHomotopy
    straightline::StraightLineHomotopy{Start, Target}
    U::Matrix{Complex128}
end

function RandomUnitaryPath(start::Systems.AbstractSystem, target::Systems.AbstractSystem)
    m, n = size(start)
    # construct a random unitary matrix
    U = qr(randn(n,n) + im * randn(n,n))[1]
    RandomUnitaryPath(Homotopies.StraightLineHomotopy(start, target), U)
end

# We have to define the size
Base.size(H::RandomUnitaryPath) = size(H.straightline)

# Cache for efficient evaluation
struct RandomUnitaryPathCache{C, T1, T2} <: Homotopies.AbstractHomotopyCache
    straightline::C
    y::Vector{T1}
    jac_buffer::Matrix{T2}
    dt_buffer::Vector{T2}
    U_t::Matrix{Complex128}
    U_t_buffer::Matrix{Complex128}
end

function Homotopies.cache(H::RandomUnitaryPath, x, t)
    y = H.U * x
    straightline = Homotopies.cache(H.straightline, y, t)
    U_t = copy(H.U)
    U_t_buffer::Matrix{Complex128} = copy(U_t)
    jac_buffer = Homotopies.jacobian(H.straightline, y, t, straightline)
    dt_buffer = jac_buffer * y
    RandomUnitaryPathCache(straightline, y, jac_buffer, dt_buffer, U_t, U_t_buffer)
end

# Evaluation and Differentation
"""
    Ut_mul_x!(cache, U, t)

Evaluate ``U(t)x`` and return the result.
"""
function Ut_mul_x!(cache, U, x, t)
    y, U_t, buffer = cache.y, cache.U_t, cache.U_t_buffer
    # We start with U * (the 2x2 sin-cos block + I)
    buffer .= U
    s, c = sin(2π*t), cos(2π*t)
    for i=1:size(U, 1)
        buffer[i, 1] = U[i,2] * s + U[i,1] * c
        buffer[i, 2] = U[i,2] * c - U[i,1] * s
    end
    A_mul_B!(cache.y, A_mul_Bc!(U_t, buffer, U), x)
end

"""
    U_dot_t_mul_x!(cache, U, t)

Evaluate ``U'(t)x`` and return the result.
"""
function U_dot_t_mul_x!(cache, U, x, t)
    y, U_t, buffer = cache.y, cache.U_t, cache.U_t_buffer
    # We start with U * (the derivative of the 2x2 sin-cos block + 0)
    buffer .= zero(eltype(U))
    s, c = 2π*sin(2π*t), 2π*cos(2π*t)
    for i=1:size(U, 1)
        buffer[i, 1] =  U[i,2] * c - U[i,1] * s
        buffer[i, 2] = -U[i,2] * s - U[i,1] * c
    end

    A_mul_B!(y, A_mul_Bc!(U_t, buffer, U), x)
end

function Homotopies.evaluate!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.evaluate!(out, H.straightline, y, t, cache.straightline)
end

function Homotopies.jacobian!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.jacobian!(cache.jac_buffer, H.straightline, y, t, cache.straightline)
    A_mul_B!(out, cache.jac_buffer, cache.U_t)
end

function Homotopies.dt!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    # chain rule
    Homotopies.jacobian_and_dt!(cache.jac_buffer, out, H.straightline, y, t, cache.straightline)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t)
    A_mul_B!(cache.dt_buffer, cache.jac_buffer, y_dot)

    out .+= cache.dt_buffer
end

function Homotopies.evaluate_and_jacobian!(val, jac, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.evaluate_and_jacobian!(val, cache.jac_buffer, H.straightline, y, t, cache.straightline)
    A_mul_B!(jac, cache.jac_buffer, cache.U_t)
end

function Homotopies.jacobian_and_dt!(jac, dt, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    # chain rule
    Homotopies.jacobian_and_dt!(cache.jac_buffer, dt, H.straightline, y, t, cache.straightline)
    A_mul_B!(jac, cache.jac_buffer, cache.U_t)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t)
    A_mul_B!(cache.dt_buffer, cache.jac_buffer, y_dot)

    dt .+= cache.dt_buffer
    nothing
end

# Test the implementation
# =======================
#
# We provide a test suite to test the homotopy interface. Let's use this
# to check our implementation.
@polyvar x y z

F = Systems.SPSystem([x^2*y-3x*z, z^2*x+3y^2])
G = Systems.SPSystem([z*x^2-3x*y^2, z^3*x-2x*y*z^2])
size(RandomUnitaryPath(G, F))
InterfaceTest.homotopy(RandomUnitaryPath(G, F))


# And now we can solve systems with our custom homotopy
solve([x^2 - y, y^3*x-x], homotopy=RandomUnitaryPath)
