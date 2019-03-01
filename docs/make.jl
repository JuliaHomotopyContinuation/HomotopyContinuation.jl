using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "Solving Polynomial Systems" => "solving.md",
        "Systems" => "systems.md",
        "Homotopies" => "homotopies.md",
        "Predictors and Correctors" => "predictors-correctors.md",
        "Path tracker" => "pathtracking.md",
        "Newton's method" => "newton.md",
        "Reference" => "reference.md"
        ],
    strict=true
)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git"
)
