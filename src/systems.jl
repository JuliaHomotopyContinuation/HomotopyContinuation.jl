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
    evaluate_and_jacobian!,
    differentiate_parameters!,
    differentiate_parameters

export NullCache

include("systems/sp_system.jl")
include("systems/fixed_homotopy.jl")
include("systems/fp_system.jl")
include("systems/totaldegree_system.jl")

# This has to be here otherwise the compiler crashes
"""
    Base.size(F::AbstractSystem)

Returns a tuple `(m, n)` indicating that `F` is a system of `m` polynomials `m` in `n` variables.
"""
Base.size(::AbstractSystem) = error("Mandatory to define `Base.size` for `AbstractSystem`s")

end
