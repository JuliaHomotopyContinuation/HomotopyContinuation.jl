module Homotopies

import ..HomotopiesBase: AbstractHomotopy,
    AbstractHomotopyCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    precondition!, update!

import ..Systems
import ..Systems: AbstractSystem, AbstractSystemCache
import ..Utilities: randomish_gamma

export AbstractHomotopy,
    AbstractHomotopyCache,
    nvariables,
    cache,
    evaluate!, evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    precondition!, update!

include("homotopies/homotopy_witch_cache.jl")

include("homotopies/straight_line.jl")
include("homotopies/fixed_point.jl")
include("homotopies/patched_homotopy.jl")
include("homotopies/patch_switcher_homotopy.jl")

end
