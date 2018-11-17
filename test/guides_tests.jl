@testset "Guides" begin
    guides_dir = joinpath(@__DIR__, "..",  "guides", "src")
    for f in readdir(guides_dir)
        if endswith(f, ".jl")
            include(joinpath(guides_dir, f))
            @test true
        end
    end
end
