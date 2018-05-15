__precompile__()

module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Endgame,
        Homotopies,
        PathTracking,
        Predictors,
        Problems,
        ProjectiveVectors,
        Solving,
        StepLength,
        Systems,
        Utilities

    include("utilities.jl")
    include("parallel.jl")
    include("projective_vectors.jl")
    include("affine_patches.jl")

    include("systems_base.jl")
    include("homotopies_base.jl")
    include("systems.jl")
    include("homotopies.jl")
    include("problems.jl")
    include("predictors.jl")
    include("correctors.jl")
    include("prediction_correction.jl")

    include("step_length.jl")

    include("path_tracking.jl")
    include("endgame.jl")

    include("solving.jl")
    include("solve.jl")


    import .Solving: nresults,
        nfinite, nsingular, natinfinity, nfailed,
        finite, results, failed, atinfinity,
        solution, residual, start_solution, issuccess,
        isfailed, isatinfinity, isfinite, issingular

    export nresults,
        nfinite, nsingular, natinfinity, nfailed,
        finite, results, failed, atinfinity,
        solution, residual, start_solution, issuccess,
        isfailed, isatinfinity, isfinite, issingular
end #
