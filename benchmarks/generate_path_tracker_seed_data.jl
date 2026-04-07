using HomotopyContinuation
using LinearAlgebra

include("../test/test_systems.jl")
include("path_tracker_benchmark_utils.jl")

using .PathTrackerBenchmarkUtils

function generate_path_tracker_seed_data(
    output_dir::AbstractString = DEFAULT_SEED_DATA_DIR;
    system::Union{Nothing,AbstractString} = nothing,
)
    output_dir = abspath(output_dir)
    systems = PathTrackerBenchmarkUtils.selected_seed_systems(
        TEST_SYSTEM_COLLECTION;
        system = system,
    )
    entries = SeedManifestEntry[]

    mkpath(output_dir)

    for (system_name, F) in systems
        println("Generating seed data for $(system_name)")
        for seed_index = 1:DEFAULT_SEED_COUNT
            seed_dir = seed_directory(output_dir, system_name, seed_index)
            mkpath(seed_dir)

            seed = monodromy_seed(system_name, seed_index)
            monodromy_result = monodromy_solve(
                F;
                seed = seed,
                duplicate_check = :certified,
                max_loops_no_progress = 5,
                target_solutions_count = 2875,
                threading = true,
                show_progress = true,
            )

            sols = solutions(monodromy_result)
            certification = certify(
                F,
                sols,
                parameters(monodromy_result);
                show_progress = false,
                threading = true,
            )
            expected_solution_count = ndistinct_certified(certification)
            found_solution_count = length(sols)
            if found_solution_count != expected_solution_count
                error(
                    "Seed $(seed_index) for $(system_name) was not fully certified: " *
                    "$(found_solution_count) certified solutions vs $(expected_solution_count) distinct certified solutions.",
                )
            end

            start_parameters_file = joinpath(seed_dir, "start_parameters.txt")
            start_solutions_file = joinpath(seed_dir, "start_solutions.txt")
            target_parameters_file = joinpath(seed_dir, "target_parameters.txt")

            write_parameters(
                start_parameters_file,
                ComplexF64.(parameters(monodromy_result)),
            )
            write_solutions(start_solutions_file, sols)
            write_parameters(
                target_parameters_file,
                sample_target_parameters(system_name, seed_index, nparameters(F)),
            )

            push!(
                entries,
                SeedManifestEntry(
                    string(system_name),
                    seed_index,
                    seed,
                    expected_solution_count,
                    nparameters(F),
                    nvariables(F),
                    relative_to_root(output_dir, start_parameters_file),
                    relative_to_root(output_dir, start_solutions_file),
                    relative_to_root(output_dir, target_parameters_file),
                ),
            )
        end
    end

    manifest_path = joinpath(output_dir, MANIFEST_FILENAME)
    write_seed_manifest(manifest_path, entries)
    println("Wrote $(length(entries)) seed entries to $(manifest_path)")
    manifest_path
end

function parse_seed_generator_args(args)
    systems = Set(first.(parameterized_test_systems(TEST_SYSTEM_COLLECTION)))
    if isempty(args)
        return (DEFAULT_SEED_DATA_DIR, nothing)
    elseif length(args) == 1
        return args[1] in systems ? (DEFAULT_SEED_DATA_DIR, args[1]) : (args[1], nothing)
    else
        return (args[1], args[2])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_dir, system = parse_seed_generator_args(ARGS)
    generate_path_tracker_seed_data(output_dir; system = system)
end
