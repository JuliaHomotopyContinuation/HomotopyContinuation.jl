struct SolverOptions
    endgame_start::Float64
    abstol::Float64
    at_infinity_tol::Float64
    singular_tol::Float64
    refinement_maxiters::Int
    verbose::Bool
    pathcrossing_tolerance::Float64
    pathcrossing_check::Bool
end

function SolverOptions(;abstol::Float64=1e-8,
    at_infinity_tol::Float64=1e-10,
    singular_tol::Float64=1e4,
    endgame_start::Float64=0.1,
    refinement_maxiters::Int=100,
    verbose=false,
    pathcrossing_tolerance=1e-8,
    pathcrossing_check=true)
    SolverOptions(
        endgame_start, abstol, refinement_maxiters,
        verbose, pathcrossing_tolerance, pathcrossing_check)
end

function is_solver_options_kwarg(kwarg)
    kwarg == :endgame_start ||
    kwarg == :abstol ||
    kwarg == :refinement_maxiters ||
    kwarg == :verbose ||
    kwarg == :pathcrossing_tolerance ||
    kwarg == :pathcrossing_check
end
