using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "Solving general systems" => "solving.md",
        "Solving paremeterized systems with monodromy" => "monodromy.md",
        "Sorting arrays of solutions" => "sorting.md",
        "PathTracker" => "path_tracker.md",
        "CoreTracker" => "core_tracker.md",
        "Newton's method" => "newton.md",
        "Homotopies" => "homotopies.md",
        "Predictors and correctors" => "predictors-correctors.md",
        "Data structures for polynomial systems" => "systems.md",
        "Norms and distances" => "norms_distances.md",
        "Reference" => "reference.md"
        ],
    strict=true)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git"
)
