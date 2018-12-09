module Homotopies

import ..HomotopiesBase: AbstractHomotopy,
    AbstractHomotopyCache,
    NullCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    basehomotopy

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import ..Utilities: randomish_gamma

export AbstractHomotopy,
    AbstractHomotopyCache,
    NullCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    basehomotopy

# This has to be here otherwise the compiler crashes
"""
    Base.size(H::AbstractHomotopy)

Returns a tuple `(m, n)` indicating that `H` is a homotopy of `m` polynomials `m` in `n` variables.
"""
Base.size(::H) where {H<:AbstractHomotopy} = error("Obligatory to define `Base.size($H)`")


include("homotopies/homotopy_witch_cache.jl")

include("homotopies/straight_line.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/fixed_point.jl")
include("homotopies/patched_homotopy.jl")
include("homotopies/patch_switcher_homotopy.jl")

end
