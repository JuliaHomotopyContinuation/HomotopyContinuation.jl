using Documenter, HomotopyContinuation2

makedocs(
    sitename = "HomotopyContinuation",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
)

deploydocs(
    repo = "github.com/saschatimme/HomotopyContinuation2.jl.git",
)
