using Documenter, HomotopyContinuation
import LinearAlgebra

makedocs(
    sitename = "HomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md",
        "Solving polynomial systems" => [
            "solving.md",
            "solver.md"
        ],
        "Solving parametrized systems with monodromy" => "monodromy.md",
        "Input" => "input.md",
        "Tracking paths" => [
            "path_tracker.md",
            "core_tracker.md",
        ],
        "Homotopies" => "homotopies.md",
        "Polynomial systems" => "systems.md",
        "Reference" => "reference.md"
        ],
    strict = false,
    doctest = false,
)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git"
)
