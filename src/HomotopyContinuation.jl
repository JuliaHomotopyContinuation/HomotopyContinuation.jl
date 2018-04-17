__precompile__()

module HomotopyContinuation

    import DynamicPolynomials: @polyvar
    export @polyvar

    export AffinePatches,
        Correctors,
        Endgame,
        NewHomotopies,
        PathTracking,
        Predictors,
        Problems,
        StepSize,
        Systems,
        Utilities

    include("utilities.jl")
    include("systems.jl")
    include("new_homotopies.jl")
    include("problems.jl")
    include("predictors.jl")
    include("correctors.jl")
    include("prediction_correction.jl")
    include("affine_patches.jl")
    include("step_length.jl")

    include("path_trackers.jl")
    include("path_tracking.jl")
    include("endgame.jl")

    include("solving.jl")
    include("solve.jl")
end #
