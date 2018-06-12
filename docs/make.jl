using Documenter, HomotopyContinuation

makedocs(
    format = :html,
    sitename = "HomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md"
        ]
)

deploydocs(
    repo   = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git",
    target = "build",
    julia = "0.6",
    osname = "linux",
    deps   = nothing,
    make   = nothing
)
