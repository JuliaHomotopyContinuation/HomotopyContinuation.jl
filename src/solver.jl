export Solver

const AH{T} = AbstractHomotopy{T}
const APA = AbstractPathtrackingAlgorithm
const AEA = AbstractEndgameAlgorithm

mutable struct Solver{
    H<:AH,
    P<:Pathtracker,
    E<:Endgamer}
    homotopy::H
    pathtracker::P
    endgamer::E

    #startvalues::StartIter

    options::SolverOptions
end

function Base.deepcopy(s::Solver)
    Solver(deepcopy(s.homotopy), deepcopy(s.pathtracker),
        deepcopy(s.endgamer), deepcopy(s.options))
end

function Solver(H::AH{T}; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H); kwargs...)
end
function Solver(H::AH{T}, pa::APA; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), pa; kwargs...)
end
function Solver(H::AH{T}, pa::APA, ea::AEA; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), pa, ea; kwargs...)
end
function Solver(H::AH{T}, pa::APA, ea::AEA, HPT::Type{<:AbstractFloat}; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), pa, ea, HPT; kwargs...)
end

function Solver(
    H::AH{Complex{T}},
    #startvalues,
    pathtracking_algorithm::PA=SphericalPredictorCorrector(),
    endgame::EA=CauchyEndgame(),
    HT=widen(T);
    apply_gammatrick=true,
    kwargs...) where {T<:AbstractFloat, PA<:APA, EA<:AEA}

    # I would love to have pathtracking_algorithm, endgame and HT as kwarg, but currently (0.6.1)
    # Julia will not dispatch on kwargs. Hence, to get type stability we need is a positional
    # argument

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
    pathtracker = Pathtracker(H, pathtracking_algorithm, HT; pathtracker_kwargs...)

    endgamer_kwargs = filter_kwargs(is_endgamer_kwarg, kwargs)
    endgamer = Endgamer(endgame, pathtracker; endgamer_kwargs...)

    Solver{typeof(H), typeof(pathtracker), typeof(endgamer)}(H, pathtracker, endgamer, solver_options)
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
