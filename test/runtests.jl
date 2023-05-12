using HomotopyContinuation, LinearAlgebra, Test, Parameters, Random
using HomotopyContinuation.DoubleDouble: ComplexDF64
import Arblib, Combinatorics
using ProjectiveVectors: PVector, affine_chart
using TreeViews: TreeViews
const HC = HomotopyContinuation
import SemialgebraicSets

set_default_compile(:none)

include("test_systems.jl")

function test_treeviews(x)
    @test TreeViews.hastreeview(x)
    @test_nowarn TreeViews.treelabel(devnull, x, MIME"application/prs.juno.inline"())
    for i = 1:TreeViews.numberofnodes(x)
        @test_nowarn TreeViews.nodelabel(devnull, x, i, MIME"application/prs.juno.inline"())
        @test_nowarn TreeViews.treenode(x, i)
    end
end

function test_show_juno(x)
    @test show(stdout, MIME("application/prs.juno.inline"), x) === x
end

hc_testset = get(ENV, "HC_TESTSET", "none")

Random.seed!(0x8b868a97)

include("homotopies_test.jl")
# include("systems_test.jl")

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
