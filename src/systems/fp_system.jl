export FPSystem

"""
    FPSystem(polynomials, vars) <: AbstractSystem

Create a polynomial system using the [`FixedPolynomials`](https://github.com/JuliaAlgebra/FixedPolynomials.jl) package.
"""
struct FPSystem{T} <: AbstractSystem
    system::FP.System{T}
end

FPSystem(polys::Vector{<:MP.AbstractPolynomial}, vars) = FPSystem(FP.System(polys, vars))
function FPSystem(polys::Vector{<:MP.AbstractPolynomial};
    variables=nothing, parameters=nothing)

    parameters === nothing || error("Parameters are not supported by FPSystem")

    if variables === nothing
        FPSystem(FP.System(polys))
    else
        FPSystem(FP.System(polys, variables))
    end
 end

struct FPSystemCache{JC<:FP.JacobianConfig} <: AbstractSystemCache
    config::JC
end

cache(F::FPSystem, x) = FPSystemCache(FP.config(F.system, x))

Base.size(F::FPSystem) = (length(F.system), FP.nvariables(F.system))

@propagate_inbounds evaluate!(u, F::FPSystem, x, c::FPSystemCache) = FP.evaluate!(u, F.system, x, c.config)
evaluate(F::FPSystem, x, c::FPSystemCache) = FP.evaluate(F.system, x, c.config)
@propagate_inbounds jacobian!(U, F::FPSystem, x, c::FPSystemCache) = FP.jacobian!(U, F.system, x, c.config)
jacobian(F::FPSystem, x, c::FPSystemCache) = FP.jacobian(F.system, x, c.config)
@propagate_inbounds function evaluate_and_jacobian!(u, U, F::FPSystem, x, c::FPSystemCache)
    FP.evaluate_and_jacobian!(u, U, F.system, x, c.config)
end
function evaluate_and_jacobian(F::FPSystem, x, c::FPSystemCache)
    FP.evaluate_and_jacobian(F.system, x, c.config)
end
