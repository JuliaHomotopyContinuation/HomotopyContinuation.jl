include("setup.jl")

@testset "HomotopyContinuation" begin
    if (hc_testset == "none") || (hc_testset == "model_kit")
        include("model_kit_test.jl")
    end
    if (hc_testset == "none") || (hc_testset == "non_model_kit")
        include("utils_test.jl")
        include("double_double_test.jl")
        include("binomial_system_test.jl")
        include("norm_test.jl")
        include("voronoi_tree_test.jl")
        include("unique_points_test.jl")
        include("linear_algebra_test.jl")
        include("systems_test.jl")
        include("homotopies_test.jl")
        include("tracker_test.jl")
        include("valuation_test.jl")
        include("linear_test.jl")
        include("endgame_tracker_test.jl")
        include("polyhedral_test.jl")
        include("solve_test.jl")
        include("endgame_test.jl")
        include("monodromy_test.jl")
        include("witness_set_test.jl")
        include("certification_test.jl")
        include("semialgebraic_sets_test.jl")
        include("nid_test.jl")
    end

    # if "extensive" in ARGS
    #     include("extensive/extensive_test.jl")
    # end
end
