using Literate, Documenter, DocumenterMarkdown

for f in readdir(joinpath(@__DIR__, "src"))
    if endswith(f, ".jl")
        Literate.markdown(joinpath(@__DIR__, "src", f), joinpath(@__DIR__, "literate_build"), credit=false)
    end
end

makedocs(format=:markdown,
    root=@__DIR__,
    source="literate_build")
