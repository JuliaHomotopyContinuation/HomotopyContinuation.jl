struct SolverOptions
    endgame_start::Float64
    abstol::Float64
    at_infinity_tol::Float64
    singular_tol::Float64
    refinement_maxiters::Int
    verbose::Bool
    pathcrossing_tolerance::Float64
    pathcrossing_check::Bool
    parallel_type::Symbol
    batch_size::Int
    apply_gammatrick::Bool
    gamma::Complex128
end

function SolverOptions(;abstol::Float64=1e-12,
    at_infinity_tol::Float64=1e-10,
    singular_tol::Float64=1e4,
    endgame_start::Float64=0.1,
    refinement_maxiters::Int=100,
    verbose=false,
    pathcrossing_tolerance=1e-8,
    pathcrossing_check=true,
    parallel_type=:pmap,
    batch_size::Int=1,
    apply_gammatrick=true,
    gamma=exp(im*rand()*2π))
    SolverOptions(
        endgame_start, abstol, at_infinity_tol, singular_tol, refinement_maxiters,
        verbose, pathcrossing_tolerance, pathcrossing_check, parallel_type, batch_size,
        apply_gammatrick, gamma)
end

function is_solver_options_kwarg(kwarg)
    kwarg == :endgame_start ||
    kwarg == :abstol ||
    kwarg == :at_infinity_tol ||
    kwarg == :singular_tol ||
    kwarg == :refinement_maxiters ||
    kwarg == :verbose ||
    kwarg == :pathcrossing_tolerance ||
    kwarg == :pathcrossing_check ||
    kwarg == :parallel_type ||
    kwarg == :batch_size ||
    kwarg == :apply_gammatrick ||
    kwarg == :gamma
end
