module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Homotopies,
        InterfaceTest,
        Predictors,
        Input,
        Problems,
        Systems,
        Utilities

    include("utilities.jl")
    include("parallel.jl")
    include("affine_patches.jl")

    include("systems_base.jl")
    include("homotopies_base.jl")
    include("systems.jl")
    include("homotopies.jl")
    include("input.jl")
    include("problems.jl")
    include("predictors.jl")
    include("correctors.jl")

    include("path_tracking.jl")
    include("endgaming.jl")

    include("solving.jl")
    include("solve.jl")
    include("monodromy.jl")

    include("interface_test.jl")

    import .Homotopies: StraightLineHomotopy, FixedPointHomotopy, ParameterHomotopy
    export StraightLineHomotopy, FixedPointHomotopy, ParameterHomotopy

    import .Systems: FPSystem, SPSystem
    export FPSystem, SPSystem

    import .Utilities: homogenize, uniquevar, ishomogenous
    export homogenize, uniquevar, ishomogenous

    import LinearAlgebra: issuccess
    export issuccess

end
