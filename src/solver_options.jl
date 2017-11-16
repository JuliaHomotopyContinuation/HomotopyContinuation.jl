struct SolverOptions
    endgame_start::Float64
    abstol::Float64
    at_infinity_tol::Float64
    singular_tol::Float64
    refinement_maxiters::Int
    verbose::Bool
end

function SolverOptions(;abstol::Float64=1e-8,
    at_infinity_tol::Float64=1e-10,
    singular_tol::Float64=1e4,
    endgame_start::Float64=0.1,
    refinement_maxiters::Int=100,
    verbose=false)
    SolverOptions(endgame_start, abstol, at_infinity_tol, singular_tol, refinement_maxiters, verbose)
end

function is_solver_options_kwarg(kwarg)
    kwarg == :endgame_start ||
    kwarg == :abstol ||
    kwarg == :refinement_maxiters ||
    kwarg == :verbose
end
