using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    format = :html,
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "Solving Polynomial Systems" => "solving.md",
        "Systems" => "systems.md",
        "Homotopies" => "homotopies.md",
        "Predictors and Correctors" => "predictors-correctors.md",
        "Path tracker" => "pathtracking.md",
        "Reference" => "reference.md"
        ]
)

deploydocs(
    repo   = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git",
    target = "build",
    julia = "0.7",
    osname = "linux",
    deps   = nothing,
    make   = nothing
)
