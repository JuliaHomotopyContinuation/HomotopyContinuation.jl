module Systems

import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticPolynomials
const SP = StaticPolynomials

import ..HomotopiesBase
import ..SystemsBase: AbstractSystem,
    AbstractSystemCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    evaluate_and_jacobian,
    evaluate_and_jacobian!

include("systems/sp_system.jl")
include("systems/fixed_homotopy.jl")

end
