module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Endgaming,
        Homotopies,
        InterfaceTest,
        Monodromy,
        PathTracking,
        Predictors,
        Input,
        Problems,
        Solving,
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

    import .Solving: AffineResult, ProjectiveResult, PathResult, solution,
        residual, startsolution, issuccess,
        isfailed, isaffine, isprojective,
        isatinfinity, issingular, isnonsingular,
        nresults, nfinite, nsingular, natinfinity, nfailed, nnonsingular,
        finite, results, mapresults, solutions, realsolutions, failed, atinfinity, singular,
        nonsingular, seed, nreal, multiplicities, uniquesolutions, statistics

    export AffineResult, ProjectiveResult, PathResult, solution,
        residual, startsolution, #issuccess,
        isfailed, isaffine, isprojective,
        isatinfinity, issingular, isnonsingular,
        nresults, nfinite, nsingular, natinfinity, nfailed, nnonsingular,
        finite, results, mapresults, solutions, realsolutions, failed, atinfinity, singular,
        nonsingular, seed, nreal, multiplicities, uniquesolutions, statistics

    import .Homotopies: StraightLineHomotopy, FixedPointHomotopy, ParameterHomotopy
    export StraightLineHomotopy, FixedPointHomotopy, ParameterHomotopy

    import .Systems: FPSystem, SPSystem
    export FPSystem, SPSystem

    import .Utilities: homogenize, uniquevar, ishomogenous
    export homogenize, uniquevar, ishomogenous

    import .PathTracking: pathtracker, pathtracker_startsolutions
    export pathtracker, pathtracker_startsolutions

    import .Monodromy: monodromy_solve
    export monodromy_solve

    import LinearAlgebra: issuccess
    export issuccess

end #
