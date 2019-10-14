using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "Homotopy Continuation",
    pages = [
        "Introduction" => "index.md",
        "Solving general systems" => "solving.md",
        "Solving parametrized systems with monodromy" => "monodromy.md",
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
