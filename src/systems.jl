module Systems

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ..HomotopiesBase
import ..SystemsBase: AbstractSystem,
    AbstractSystemCache,
    NullCache,
    cache,
    evaluate,
    evaluate!,
    jacobian,
    jacobian!,
    evaluate_and_jacobian,
    evaluate_and_jacobian!

export NullCache

include("systems/sp_system.jl")
include("systems/fixed_homotopy.jl")
include("systems/fp_system.jl")

end
