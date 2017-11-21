export Solver

const AH{T} = AbstractHomotopy{T}
const APA = AbstractPathtrackingAlgorithm
const AEA = AbstractEndgameAlgorithm

"""
    Solver(homotopy, pathtracking_algorithm=SphericalPredictorCorrector(), endgame=CauchyEndgame(); kwargs...)

Create a mutable `Solver` struct. This contains a
[`Pathtracker`](@ref) and an [`Endgamer`](@ref), everything you need to solve the given
homotopy. `Solver` supports the following options:

* `endgame_start=0.1`: Where the endgame starts
* `abstol=1e-12`: The desired accuracy of the final roots
* `at_infinity_tol=1e-10`: An point is at infinity if the maginitude of the homogenous variable
is less than `at_infinity_tol`.
* `singular_tol=1e4`: If the winding number is 1 but the condition number is larger than
`singular_tol` then the root is declared as singular.
* `refinement_maxiters=100`: The maximal number of newton iterations to achieve `abstol`.
* `verbose=false`: Print additional warnings / informations
* `apply_gammatrick=true`: This modifies the start system to make it generic.
* `gamma=apply_gammatrick ? exp(im*2π*rand()) : complex(1.0)`: You can overwrite the default gamma.
    This is useful if you want to rerun only some paths.
* `pathcrossing_tolerance=1e-8`: The tolerance for when two paths are considered to be crossed.
* `pathcrossing_check=true`: Enable the pathcrossing check.
* `parallel_type=:pmap`: Currently there are two modes: `:pmap` will use `pmap` for parallelism
and `:none` will use the standard `map`. `:pmap` is by defautl enabled since it works reliable,
but if you develop new algorithms you probably want to disable parallelism.
* `batch_size=1`: The `batch_size` for `pmap` if `parallel_type` is `:pmap`.

For instance, to solve the homotopy `H` with starting values `s` with no endgame and a singular tolerance of 1e5, write

```julia
    solve(H, s, endgame_start=0.0, singular_tol=1e5)
```

To solve the polynomial system ``f`` with the same options write

```julia
    solve(f, endgame_start=0.0, singular_tol=1e5)
```
"""
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

# Base.copy(s::Solver) = deepcopy(s)
# function Base.deepcopy(s::Solver)
#     Solver(deepcopy(s.homotopy), deepcopy(s.pathtracker),
#         deepcopy(s.endgamer),  deepcopy(s.options))
# end

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

    solver_options_kwargs = filter_kwargs(is_solver_options_kwarg, kwargs)
    solver_options = SolverOptions(;solver_options_kwargs...)

    γ = one(Complex128)
    if solver_options.apply_gammatrick
        γ = solver_options.gamma
    end
    H = convert(typeof(H), gammatrick(H, convert(Complex{T}, γ)))


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
