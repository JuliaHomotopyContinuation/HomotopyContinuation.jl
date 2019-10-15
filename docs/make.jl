using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "The solve function" => "solving.md",
        "Solving parametrized systems with monodromy" => "monodromy.md",
        "The solver struct" => "solver.md",
        "PathTracker" => "path_tracker.md",
        "CoreTracker" => "core_tracker.md",
        "Homotopies" => "homotopies.md",
        "Polynomial systems" => "systems.md",
        "Reference" => "reference.md"
        ],
    strict=false)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git"
)
