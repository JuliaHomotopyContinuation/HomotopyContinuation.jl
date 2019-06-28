@testset "Polyhedral" begin
    @testset "Start Solutions Iterator" begin
        f = equations(cyclic(5))
        iter = HC.PolyhedralStartSolutionsIterator(f)

        mv = 0
        for (cell, X) in iter
            mv += cell.volume
            @test size(X,2) == cell.volume
            @test all(1:size(X,2)) do j
                x = @view X[:,j]
                all(1:length(cell.indices)) do i
                    a, b = cell.indices[i]
                    v_i = iter.start_coefficients[i][a] * prod(x.^iter.support[i][:,a]) +
                          iter.start_coefficients[i][b] * prod(x.^iter.support[i][:,b])
                    abs(v_i) < 1e-12
                end
            end
        end
        @test mv == 70

        Random.seed!(130793) # this seed should result in the use of 128 bit hnf
        f = equations(PolynomialTestSystems.tritangents())
        iter = HC.PolyhedralStartSolutionsIterator(f)
        @test sum(cell_X -> first(cell_X).volume, iter) == 12636
    end

    @testset "Polyhedral solve" begin
        f = equations(cyclic(5))
        result = solve(f; start_system=:polyhedral)
        @test nfinite(result) == 70

        @polyvar x y
        result1 = solve([(x-3),(y-2)], start_system=:polyhedral, system_scaling=nothing)
        @test [3,2] ≈ solutions(result1)[1] atol=1e-8
        result1 = solve([(x-3),(y-2)], start_system=:polyhedral, system_scaling=:equations_and_variables)
        @test [3,2] ≈ solutions(result1)[1] atol=1e-8
    end

    @testset "All affine solutions" begin
        @polyvar x y
        f = [2y + 3y^2 - x*y^3, x + 4*x^2 - 2*x^3*y]
        res = solve(f, start_system=:polyhedral, only_torus=true)
        @test nsolutions(res) == 3
        @test ntracked(res) == 3

        res = solve(f, start_system=:polyhedral, affine_tracking=false, only_torus=false)
        @test nsolutions(res) == 6
        @test ntracked(res) == 8

        res = solve(f, start_system=:polyhedral, affine_tracking=true, only_torus=false)
        @test nsolutions(res) == 6
        @test ntracked(res) == 8
    end

    @testset "Overflow error message" begin
        f = equations(cyclooctane())
        F = [f; randn(2, 18) * [variables(f);1]]
        @test_throws OverflowError solve(F; start_system=:polyhedral)
    end
end
