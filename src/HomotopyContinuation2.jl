module HomotopyContinuation2

export ModelKit

import LinearAlgebra
import StructArrays

const LA = LinearAlgebra

include("utils.jl")

include("ModelKit.jl")
import .ModelKit

include("linear_algebra.jl")

end
