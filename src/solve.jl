export solve

"""
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.


    solve(f::Vector{<:MP.AbstractPolynomial{T}}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)
    solve(f::MP.AbstractPolynomial{T}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)

Solves the polynomial system `f` via homotopy continuation. This uses a totaldegree homotopy
of the given `homotopytype` and uses the given `algorithm` for pathtracking.
"""
function solve end

function solve(
    H::AbstractHomotopy{Complex{T}},
    startvalues;
    apply_gammatrick=true,
    pathtracking_algorithm=SphericalPredictorCorrector(),
    endgame=CauchyEndgame(),
    endgame_start::Float64=0.1,
    abstol::Float64=1e-8,
    refinement_maxiters::Int=100,
    verbose::Bool=false,
    high_precision_type = _widen(T),
    maxiters::Int=10_000,
    consecutive_successfull_steps_until_steplength_increase::Int=3,
    steplength_increase_factor::Float64=2.0,
    steplength_decrease_factor::Float64=inv(steplength_increase_factor),
    path_precision::Float64=1e-6,
    corrector_maxiters::Int=3,
    initial_steplength::Float64=0.1,
    geometric_series_factor=0.5,
    max_winding_number=8) where {T}

    if apply_gammatrick
        gammatrick!(H)
    end

    pathtracker = initialize(pathtracking_algorithm, H,
        highprecisiontype = high_precision_type,
        maxiters = maxiters,
        verbose = verbose,
        consecutive_successfull_steps_until_steplength_increase = consecutive_successfull_steps_until_steplength_increase,
        steplength_increase_factor = steplength_increase_factor,
        steplength_decrease_factor = steplength_decrease_factor,
        abstol = path_precision,
        corrector_maxiters = corrector_maxiters,
        initial_steplength = initial_steplength)

    endgamer = initialize(
        endgame,
        pathtracker,
        geometric_series_factor=geometric_series_factor,
        max_winding_number=max_winding_number)

    options = SolverOptions(endgame_start, abstol, refinement_maxiters, verbose)

    solver = Solver(H, pathtracker, endgamer, startvalues, options)
    solve(solver)
end

function solve(
    f::MP.AbstractPolynomial{T};
    homotopytype=StraightLineHomotopy,
    kwargs...) where {T<:Number}
      H, s = totaldegree(homotopytype, [f])
      solve(H, s; kwargs...)
end
function solve(f::Vector{<:MP.AbstractPolynomial{T}};
    homotopytype=StraightLineHomotopy,
    kwargs...) where {T<:Number}
     H, s = totaldegree(homotopytype, f)
     solve(H, s; kwargs...)
end
