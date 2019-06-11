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
    end

    @testset "Polyhedral solve" begin
        f = equations(cyclic(5))
        result = solve(f; start_system=:polyhedral)
        nfinite(result) == 70
    end
end
