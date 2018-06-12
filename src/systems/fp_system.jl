import FixedPolynomials
const FP = FixedPolynomials

export FPSystem, differentiate

"""
    FPSystem(polynomials, vars) <: AbstractSystem

Create a system using the `FixedPolynomials` package.
"""
struct FPSystem{T} <: AbstractSystem
    system::FP.System{T}
end

FPSystem(polys::Vector{<:MP.AbstractPolynomial}, vars) = FPSystem(FP.System(polys, vars))
FPSystem(polys::Vector{<:MP.AbstractPolynomial}) = FPSystem(FP.System(polys))

struct FPSystemCache{JC<:FP.JacobianConfig} <: AbstractSystemCache
    config::JC
end

cache(F::FPSystem, x) = FPSystemCache(FP.config(F.system, x))

Base.size(F::FPSystem) = (length(F.system), FP.nvariables(F.system))
Base.length(F::FPSystem) = length(F.system)

evaluate!(u, F::FPSystem, x, c::FPSystemCache) = FP.evaluate!(u, F.system, x, c.config)
evaluate(F::FPSystem, x, c::FPSystemCache) = FP.evaluate(F.system, x, c.config)
jacobian!(U, F::FPSystem, x, c::FPSystemCache) = FP.jacobian!(U, F.system, x, c.config)
jacobian(F::FPSystem, x, c::FPSystemCache) = FP.jacobian(F.system, x, c.config)
function evaluate_and_jacobian!(u, U, F::FPSystem, x, c::FPSystemCache)
    FP.evaluate_and_jacobian!(u, U, F.system, x, c.config)
end
function evaluate_and_jacobian(F::FPSystem, x, c::FPSystemCache)
    FP.evaluate_and_jacobian(F.system, x, c.config)
end
