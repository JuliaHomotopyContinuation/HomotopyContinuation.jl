export Solver

const AH{T} = AbstractHomotopy{T}
const APA = AbstractPathtrackingAlgorithm
const AEA = AbstractEndgameAlgorithm

mutable struct Solver{
    H<:AH,
    P<:Pathtracker,
    E<:Endgamer,
    StartIter}
    homotopy::H
    pathtracker::P
    endgamer::E

    startvalues::StartIter

    options::SolverOptions
end

function Solver(
    H::AH{Complex{T}},
    startvalues,
    pathtracking_algorithm::PA=SphericalPredictorCorrector(),
    endgame::EA=CauchyEndgame(),
    HT=widen(T);
    apply_gammatrick=true,
    kwargs...) where {T<:AbstractFloat, PA<:APA, EA<:AEA}

    # I love to have pathtracking_algorithm, endgame and HT as kwarg, but currently (0.6.1)
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
    pathtracker = Pathtracker(
        pathtracking_algorithm, H, first(startvalues), 1.0,
        solver_options.endgame_start, HT; pathtracker_kwargs...)

    endgamer_kwargs = filter_kwargs(is_endgamer_kwarg, kwargs)
    endgamer = Endgamer(endgame, pathtracker; endgamer_kwargs...)

    Solver{typeof(H), typeof(pathtracker), typeof(endgamer), typeof(startvalues)}(H, pathtracker, endgamer, startvalues, solver_options)
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


# Convenience constructors
Solver(f::MP.AbstractPolynomial; kwargs...) = Solver([f]; kwargs...)
Solver(f::MP.AbstractPolynomial, pa::APA; kwargs...) = Solver([f], pa; kwargs...)
Solver(f::MP.AbstractPolynomial, pa::APA, ea::AEA; kwargs...) = Solver([f], pa, ea; kwargs...)
Solver(f::MP.AbstractPolynomial, HT; kwargs...) = Solver([f], HT; kwargs...)
Solver(f::MP.AbstractPolynomial, HT, pa::APA; kwargs...) = Solver([f], HT, pa; kwargs...)
Solver(f::MP.AbstractPolynomial, HT, pa::APA, ea::AEA; kwargs...) = Solver([f], HT, pa, ea; kwargs...)

function Solver(F::Vector{<:MP.AbstractPolynomial{T}}; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      Solver(H, s; kwargs...)
end
function Solver(F::Vector{<:MP.AbstractPolynomial{T}}, pa::APA; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      Solver(H, s, pa; kwargs...)
end
function Solver(F::Vector{<:MP.AbstractPolynomial{T}}, pa::APA, ea::AEA; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      Solver(H, s, pa, ea; kwargs...)
end
function Solver(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      Solver(H, s; kwargs...)
end
function Solver(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}, pa::APA; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      Solver(H, s, pa; kwargs...)
end
function Solver(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}, pa::APA, ea::AEA; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      Solver(H, s, pa, ea; kwargs...)
end

function Solver(H::AH{T}, s; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), s; kwargs...)
end
function Solver(H::AH{T}, s, pa::APA; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), s, pa; kwargs...)
end
function Solver(H::AH{T}, s, pa::APA, ea::AEA; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), s, pa, ea; kwargs...)
end
function Solver(H::AH{T}, s, pa::APA, ea::AEA, HPT::Type{<:AbstractFloat}; kwargs...) where {T<:Union{Complex{<:Integer}, Real}}
    HT = promote_type(typeof(H), promote_type(Complex128, T))
    Solver(convert(HT, H), s, pa, ea, HPT; kwargs...)
end
