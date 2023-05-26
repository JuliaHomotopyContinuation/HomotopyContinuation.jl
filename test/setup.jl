using HomotopyContinuation, LinearAlgebra, Test, Parameters, Random
using HomotopyContinuation.DoubleDouble: ComplexDF64
import Arblib, Combinatorics
using ProjectiveVectors: PVector, affine_chart
using TreeViews: TreeViews
const HC = HomotopyContinuation
import SemialgebraicSets
import MixedSubdivisions

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
