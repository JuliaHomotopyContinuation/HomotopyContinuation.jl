module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Endgaming,
        Homotopies,
        InterfaceTest,
        PathTracking,
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

    import .PathTracking: PathTracker, PathTrackerResult, pathtracker, pathtracker_startsolutions,
            allowed_keywords, track, track!, setup!, iterator!,
            currx, currt, currΔt, curriters, currstatus, tol, corrector_maxiters,
            refinement_tol, refinement_maxiters, set_tol!,
            set_corrector_maxiters!, set_refinement_tol!, set_refinement_maxiters!
    export PathTracker, PathTrackerResult, pathtracker, pathtracker_startsolutions,
            allowed_keywords, track, track!, setup!, iterator!,
            currx, currt, currΔt, curriters, currstatus, tol, corrector_maxiters,
            refinement_tol, refinement_maxiters, set_tol!,
            set_corrector_maxiters!, set_refinement_tol!, set_refinement_maxiters!

    import LinearAlgebra: issuccess
    export issuccess

end
