import MultivariatePolynomials
const MP=MultivariatePolynomials
export ParameterHomotopy, ParameterHomotopyCache

"""
    ParameterHomotopy(f, variables, parameters, start, target; gamma=exp(i * 2π*rand()))

Construct the homotopy ``H(x, t) = F(x, t * start + (1-t) * target)``, where `start` and `target` are a vector of parameters of ``F``. The input ``parameters`` specifies the parameter variables of ``F``. Neccessarily, ``length(parameters) == length(start) == length(target)``.
"""
struct ParameterHomotopy{T<:AbstractSystem, S1<:Number, S2<:Number} <: AbstractHomotopy
    f::T
    npolys::Int

    variables::Vector{Int}
    parameters::Vector{Int}
    nvariables::Int
    nparameters::Int

    start::Vector{S1}
    target::Vector{S2}
    gamma::Complex{Float64}
end
function ParameterHomotopy(
    f::AbstractSystem,
    variables::Vector{Int},
    parameters::Vector{Int},
    start::Vector{S1},
    target::Vector{S2}
    ; gamma=randomish_gamma()) where {S1<:Number, S2<:Number}


    ParameterHomotopy(f, length(f), variables, parameters, length(variables), length(parameters), start, target, gamma)
end
(H::ParameterHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

"""
    ParameterHomotopyCache

An simple cache for `ParameterHomotopy`s consisting of the caches for the
start and target system as well as a `Vector` and a `Matrix`.
"""
struct ParameterHomotopyCache{fC, T, S} <: AbstractHomotopyCache
    f::fC
    z::Vector{T} # z is the vector z = [x; ta + (1-t)b]
    U::Matrix{S} # for storing the jacobian
    U_t::Union{Vector{S}, Matrix{S}} # for storing the partial derivatives wrt t
end

function cache(H::ParameterHomotopy, x, t)

    z = [x; t * H.start + (1-t) * H.target]
    f_cache = Systems.cache(H.f, z)
    U = Systems.jacobian(H.f, z, f_cache)
    U_t = U[:, H.parameters]

    ParameterHomotopyCache(f_cache, z, U, U_t)
end

Base.size(H::ParameterHomotopy) = (H.npolys, H.nvariables)

"""
    gamma(H::ParameterHomotopy)

Obtain the gamma used in the ParameterHomotopy.
"""
Base.Math.gamma(H::ParameterHomotopy) = H.gamma

"""
    γ(H)

Obtain the gamma used in the ParameterHomotopy.
"""
γ(H::ParameterHomotopy) = gamma(H)

function update_z!(c::ParameterHomotopyCache, H, x, t)
    for i in H.variables
        c.z[i] = x[i]
    end
    for (i,j) in enumerate(H.parameters)
        c.z[j] = t .* H.start[i] + (1-t) .* H.target[i]
    end
end

function evaluate!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    update_z!(c, H, x, t)
    Systems.evaluate!(u, H.f, c.z, c.f)

    u
end

function dt!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    update_z!(c, H, x, t)
    Systems.jacobian!(c.U, H.f, c.z, c.f)
    for i in 1:H.npolys
        for (j1, j2) in enumerate(H.parameters)
            c.U_t[i,j1] = c.U[i,j2]
        end
    end
    A_mul_B!(u, c.U_t, H.start - H.target)

    u
end

function jacobian!(U, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    update_z!(c, H, x, t)
    Systems.jacobian!(c.U, H.f, c.z, c.f)
    for i in 1:H.npolys
        for j in H.variables
            U[i,j] = c.U[i,j]
        end
    end

    U
end

function evaluate_and_jacobian!(u, U, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    update_z!(c, H, x, t)
    Systems.evaluate_and_jacobian!(u, c.U, H.f, c.z, c.f)
    for i in 1:H.npolys
        for j in H.parameters
            U[i,j] = c.U[i,j]
        end
    end

    nothing
end

function jacobian_and_dt!(U, u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    update_z!(c, H, x, t)

    Systems.evaluate_and_jacobian!(u, c.U, H.f, c.z, c.f)
    for i in 1:H.npolys
        for j in H.variables
            U[i,j] = c.U[i,j]
        end
    end
    for i in 1:H.npolys
        for j in H.parameters
            c.U_t[i,j] = c.U[i,j]
        end
    end
    A_mul_B!(u, c.U_t, H.start - H.target)

    nothing
end

function evaluate(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    evaluate!(similar(c.U_t, H.npolys), H, x, t, c)
end
function jacobian(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    jacobian!(similar(c.U_t, H.npolys, H.nvariables), H, x, t, c)
end
function dt(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    dt!(similar(c.U_t, H.npolys), H, x, t, c)
end
