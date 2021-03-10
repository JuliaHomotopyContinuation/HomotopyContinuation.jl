using HomotopyContinuation, LinearAlgebra, Test, Parameters, Random
import Arblib
using StaticArrays
using ProjectiveVectors: PVector, affine_chart
using TreeViews: TreeViews
const HC = HomotopyContinuation
import SemialgebraicSets

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

@testset "HomotopyContinuation" begin
    include("utils_test.jl")
    include("double_double_test.jl")
    include("binomial_system_test.jl")
    include("norm_test.jl")
    include("voronoi_tree_test.jl")
    include("unique_points_test.jl")
    include("model_kit_test.jl")
    include("linear_algebra_test.jl")
    include("homotopies_test.jl")
    include("systems_test.jl")
    include("tracker_test.jl")
    include("valuation_test.jl")
    include("endgame_tracker_test.jl")
    include("linear_test.jl")
    include("polyhedral_test.jl")
    include("solve_test.jl")
    include("endgame_test.jl")
    include("monodromy_test.jl")
    include("witness_set_test.jl")
    include("certification_test.jl")
    include("semialgebraic_sets_test.jl")

    if "extensive" in ARGS
        include("extensive/extensive_test.jl")
    end
end
