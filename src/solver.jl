export Solver

mutable struct Solver{
    H<:AbstractHomotopy,
    P<:Pathtracker,
    E<:Endgamer,
    StartIter}
    homotopy::H
    pathtracker::P
    endgamer::E

    startvalues::StartIter

    options::SolverOptions
end


function Solver(H::AbstractHomotopy{Complex{T}}, startvalues, HT=widen(T);
    pathtracking_algorithm=SphericalPredictorCorrector(),
    endgame=CauchyEndgame(), kwargs...) where T

    Solver(H, startvalues, pathtracking_algorithm, endgame, HT; kwargs...)
end


function Solver(
    H::AbstractHomotopy{Complex{T}},
    startvalues,
    pathtracking_algorithm::AbstractPathtrackingAlgorithm,
    endgame::AbstractEndgameAlgorithm,
    ::Type{HT};
    apply_gammatrick=true,
    kwargs...) where {T, HT}

    # Due to some internal limitation of Julia (as of 0.6.1) we cannot have more than
    # 11 kwargs without losing proper type inference. As a workaround we can pass
    # the kwargs through to the different constructors. One problem with this approach
    # is that invalid kwargs will not throw an error (due to the filter functions.).
    # To fix this behaviour, we manually check for any invalid kwargs
    assert_valid_kwargs(kwargs)

    if apply_gammatrick
        gammatrick!(H)
    end

    solver_options_kwargs = filter_kwargs(is_solver_options_kwarg, kwargs)
    solver_options = SolverOptions(;solver_options_kwargs...)

    pathtracker_kwargs = filter_kwargs(is_pathtracker_kwarg, kwargs)
    pathtracker = Pathtracker(
        pathtracking_algorithm, H, first(startvalues), 1.0,
        solver_options.endgame_start, HT; pathtracker_kwargs...)

    endgamer_kwargs = filter_kwargs(is_endgamer_kwarg, kwargs)
    endgamer = Endgamer(endgame, pathtracker; endgamer_kwargs...)

    Solver(H, pathtracker, endgamer, startvalues, solver_options)
end

function assert_valid_kwargs(kwargs)
    for kwarg in kwargs
        kw = first(kwarg)
        is_valid =
            is_solver_options_kwarg(kw) ||
            is_pathtracker_kwarg(kw) ||
            is_endgamer_kwarg(kw)

        if !is_valid
            throw(ArgumentError("Unknown keyword argument `$(kw)` passed to `solve`"))
        end
    end
    nothing
end
