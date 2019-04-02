using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "Solving Polynomial Systems" => "solving.md",
        "Solving Systems with Monodromy" => "monodromy.md",
        "Newton's method" => "newton.md",
        "Sorting arrays of solutions" => "sorting.md",
        "Norms and Distances" => "norms_distances.md",
        "Data structures for polynomial systems" => "systems.md",
        "Homotopies" => "homotopies.md",
        "Predictors and Correctors" => "predictors-correctors.md",
        "Path tracker" => "pathtracking.md",
        "Reference" => "reference.md"
        ],
    strict=true
)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git"
)
