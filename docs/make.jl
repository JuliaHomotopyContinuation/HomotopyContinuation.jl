using Documenter, HomotopyContinuation2

makedocs(
    sitename = "HomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md",
        "Problem formulation" => [
            "ModelKit" => "model_kit.md",
            "Linear and Affine Subspaces" => "linear_affine.md",
            "Systems and Homotopies" => "systems_homotopies.md",
        ],
        "Trackers" => ["PathTracker" => "path_tracker.md", "Trackers" => "tracker.md"],
        "Miscellaneous" => "misc.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
)

deploydocs(
    repo = "github.com/saschatimme/HomotopyContinuation2.jl.git",
    push_preview = false,
)
