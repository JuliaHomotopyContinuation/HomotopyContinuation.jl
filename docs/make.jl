using Documenter, HomotopyContinuation

makedocs(
    sitename = "HomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md",
        "Problem formulation" => [
            "ModelKit" => "model_kit.md",
            "Systems" => "systems.md",
            "Homotopies" => "homotopies.md",
            "Linear Subspaces" => "linear_subspaces.md",
        ],
        "Solving Systems" => [
            "Solve (finitely many solutions)" => "solve.md",
            "Positive dimensional solution sets" => "witness_sets.md",
            "Monodromy" => "monodromy.md",
            "Results" => "result.md",
            "Certification" => "certification.md",
            "Examples" => "solve_examples.md",
            "Start systems" => "start_systems.md",
        ],
        "Trackers" => ["Tracker" => "tracker.md", "EndgameTracker" => "endgame_tracker.md"],
        "Miscellaneous" => "misc.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    strict = false,
)

deploydocs(
    repo = "github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl.git",
    push_preview = false,
)
