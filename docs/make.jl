using Documenter, HomotopyContinuation

makedocs(
    sitename = "HomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md",
        "Problem formulation" => [
            "ModelKit" => "model_kit.md",
            "Systems and Homotopies" => "systems_homotopies.md",
            "Linear and Affine Subspaces" => "linear_affine.md",
        ],
        "Solving systems" => [
            "Solve" => "solve.md",
            "Monodromy" => "monodromy.md",
            "Results" => "result.md",
            "Examples" => "solve_examples.md",
            "Start systems" => "start_systems.md",
        ],
        "Trackers" => ["Tracker" => "tracker.md", "EndgameTracker" => "endgame_tracker.md"],
        "Miscellaneous" => "misc.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git",
    push_preview = false,
)
