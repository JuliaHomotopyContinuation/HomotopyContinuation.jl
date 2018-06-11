import StaticPolynomials
const SP = StaticPolynomials

export SPSystem

"""
    SPSystem(polynomials, vars) <: AbstractSystem

Create a system using the `StaticPolynomials` package.
"""
struct SPSystem{S<:SP.AbstractSystem} <: AbstractSystem
    system::S
end

SPSystem(polys::Vector{<:MP.AbstractPolynomial}, vars) = SPSystem(SP.system(polys, vars))
SPSystem(polys::Vector{<:MP.AbstractPolynomial}) = SPSystem(SP.system(polys))

Base.size(F::SPSystem) = (SP.npolynomials(F.system), SP.nvariables(F.system))
Base.length(F::SPSystem) = SP.npolynomials(F.system)

evaluate!(u, F::SPSystem, x) = SP.evaluate!(u, F.system, x)
evaluate(F::SPSystem, x) = SP.evaluate(F.system, x)
jacobian!(U, F::SPSystem, x) = SP.jacobian!(U, F.system, x)
jacobian(F::SPSystem, x) = SP.jacobian(F.system, x)
evaluate_and_jacobian!(u, U, F::SPSystem, x) = SP.evaluate_and_jacobian!(u, U, F.system, x)
evaluate_and_jacobian(F::SPSystem, x) = SP.evaluate_and_jacobian(F.system, x)
